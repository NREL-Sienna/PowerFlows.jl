struct HomotopyHessian
    # PERF: data is stored in triplicate: here, inside pfResidual, and inside J.
    data::ACPowerFlowData
    pfResidual::ACPowerFlowResidual
    J::ACPowerFlowJacobian
    PQ_V_mags::BitVector # true iff that coordinate in the state vector is V_mag at a PQ bus
    grad::Vector{Float64}
    Hv::SparseMatrixCSC{Float64, Int32}
    t_k_ref::Base.RefValue{Float64}
end

"""Does `A += B' * B`, in a way that preserves the sparse structure of `A`, if possible.
A workaround for the fact that Julia seems to run `dropzeros!(A)` automatically if I just 
do `A .+= B' * B`."""
function A_plus_eq_BT_B!(A::SparseMatrixCSC, B::SparseMatrixCSC)
    M = B' * B
    @assert M.colptr == A.colptr && M.rowval == A.rowval
    A.nzval .+= M.nzval
    return
end

"""Compute value of gradient and Hessian at x."""
function (hess::HomotopyHessian)(x::Vector{Float64}, time_step::Int)
    t_k = hess.t_k_ref[]
    hess.pfResidual(x, time_step)
    Rv = hess.pfResidual.Rv
    hess.J(time_step)
    Jv = hess.J.Jv
    old_row, old_col = copy(hess.Hv.rowval), copy(hess.Hv.colptr)
    _update_hessian_matrix_values!(hess.Hv, Rv, hess.data, time_step)
    A_plus_eq_BT_B!(hess.Hv, Jv)
    SparseArrays.nonzeros(hess.Hv) .*= t_k
    for (bus_ix, bt) in enumerate(get_bus_type(hess.data)[:, time_step])
        if bt == PSY.ACBusTypes.PQ
            hess.Hv[2 * bus_ix - 1, 2 * bus_ix - 1] += (1 - t_k)
        end
    end
    hess.grad .= (1 - t_k) * hess.PQ_V_mags .* (x - ones(size(x, 1))) + t_k * Jv' * Rv
    @assert hess.Hv.rowval == old_row && hess.Hv.colptr == old_col
    return
end

function F_value(hess::HomotopyHessian, x::Vector{Float64}, time_step::Int)
    t_k = hess.t_k_ref[]
    hess.pfResidual(x, time_step)
    Rv = hess.pfResidual.Rv
    φ_vector = x[hess.PQ_V_mags] .- 1.0
    F_value = (1 - t_k) * 0.5 * dot(φ_vector, φ_vector) + t_k * 0.5 * dot(Rv, Rv)
    return F_value
end

function gradient_value(hess::HomotopyHessian, x::Vector{Float64}, time_step::Int)
    t_k = hess.t_k_ref[]
    hess.pfResidual(x, time_step)
    hess.J(time_step)
    Jv = hess.J.Jv
    mask = hess.PQ_V_mags
    grad = (1 - t_k) * (mask .* (x - ones(size(x, 1)))) + t_k * Jv' * hess.pfResidual.Rv
    return grad
end

function homotopy_x0(data::ACPowerFlowData, time_step::Int)
    x = calculate_x0(data, time_step)
    for (bus_ix, bt) in enumerate(get_bus_type(data)[:, time_step])
        if bt == PSY.ACBusTypes.PQ
            x[2 * bus_ix - 1] = 1.0
        end
    end
    return x
end

function HomotopyHessian(data::ACPowerFlowData, time_step::Int)
    pfResidual = ACPowerFlowResidual(data, time_step)
    Hv = _create_hessian_matrix_structure(data, time_step)
    J = ACPowerFlowJacobian(data, time_step)
    nbuses = size(get_bus_type(data), 1)
    PQ_mask = get_bus_type(data)[:, time_step] .== (PSY.ACBusTypes.PQ,)
    PQ_V_mags = collect(Iterators.flatten(zip(PQ_mask, falses(nbuses))))
    return HomotopyHessian(data, pfResidual, J, PQ_V_mags, zeros(2 * nbuses), Hv, Ref(0.0))
end

function _create_hessian_matrix_structure(data::ACPowerFlowData, time_step::Int64)
    rows = Int32[]      # I
    columns = Int32[]   # J
    values = Float64[]  # V

    # an over-estimate: I want ordered pairs of vertices that are 2 or fewer
    # steps apart, whereas this counts directed paths of 2 edges.
    numEdgePairs = sum(x -> length(x)^2, get_branch_lookup(data))
    sizehint!(rows, 4 * numEdgePairs)
    sizehint!(columns, 4 * numEdgePairs)
    sizehint!(values, 4 * numEdgePairs)

    # H is J'*J + c*I
    # J' * J is dot products of pairs of columns.
    # so look at pairs of columns and check if there's a row in which both are nonzero.
    # i.e. look at pairs of buses and see if they have a neighbor in common.

    enum_nbhrs = enumerate(data.neighbors)
    for (b1, b2) in Iterators.product(enum_nbhrs, enum_nbhrs)
        if !isdisjoint(b1[2], b2[2])
            bus_from, bus_to = b1[1], b2[1]
            # PERF: J.Jv^T * J.Jv would have fewer nonzero entries if J.Jv's sparse
            # structure took into account the bus type
            for (i, j) in Iterators.product((2 * bus_from - 1, 2 * bus_from),
                (2 * bus_to - 1, 2 * bus_to))
                push!(rows, i)
                push!(columns, j)
                push!(values, 0.0)
            end
        end
    end
    return SparseArrays.sparse(rows, columns, values)
end

"""Sets Hv equal to `F_1(x) H_{F_1}(x) + ...+ F_{2n}(x) H_{F_{2n}}(x)`,
where F_k denotes the kth power balance equation and `H_{F_k}` its Hessian.
This isn't the full Hessian of our function: it's only the terms in that come
from the second derivatives of our power balance equations. (There's also a `J'*J` term.)

What's the sparse structure of that expression? It's split into 2x2 blocks, each 
corresponding to a pair of buses. The sparse structure of a block for a pair of buses 
connected by a branch is:
   | REF| PV | PQ 
---+----+----+----
REF|    |    |   
   |    |    |
---+----+----+----
PV |    |    |
   |    |   .| . .
---+----+----+----
PQ |    |   .| . .
   |    |   .| . .
Diagonal blocks follow the same pattern as above (as if each bus is its own neighbor).
Off-diagonal blocks for a pair of buses not connected by a branch are structurally zero.
"""
function _update_hessian_matrix_values!(
    Hv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    F_value::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64)
    Yb = data.power_network_matrix.data
    Vm = view(data.bus_magnitude, :, time_step)
    θ = view(data.bus_angles, :, time_step)
    num_buses = first(size(data.bus_type))
    SparseArrays.nonzeros(Hv) .= 0.0
    for i in 1:num_buses
        bt_i = data.bus_type[i, time_step]
        bus_neighbors = data.neighbors[i]
        for k in data.neighbors[i]
            if i == k
                # ∂²Δ{Pᵢ, Qᵢ}/∂²θᵢ: PQ and PV
                if (bt_i == PSY.ACBusTypes.PQ) || (bt_i == PSY.ACBusTypes.PV)
                    Pi_θiθi =
                        Vm[i] * sum( # = -Vm[i] * Qi_Viθi
                            Vm[l] * (
                                -real(Yb[i, l]) * cos(θ[i] - θ[l])
                                -
                                imag(Yb[i, l]) * sin(θ[i] - θ[l])
                            ) for l in bus_neighbors if l != i
                        )
                    Qi_θiθi =
                        Vm[i] * sum( # = Vm[i] * Pi_Viθi
                            Vm[l] * (
                                -real(Yb[i, l]) * sin(θ[i] - θ[l])
                                +
                                imag(Yb[i, l]) * cos(θ[i] - θ[l])
                            ) for l in bus_neighbors if l != i
                        )
                    θiθis = Pi_θiθi * F_value[2 * i - 1] + Qi_θiθi * F_value[2 * i]
                    Hv[2 * i, 2 * i] += θiθis
                end
                # ∂²Δ{Pᵢ, Qᵢ}/∂Vᵢ∂θᵢ and ∂²Δ{Pᵢ, Qᵢ}/∂²Vᵢ: PQ only.
                if bt_i == PSY.ACBusTypes.PQ
                    Pi_Viθi = sum(
                        Vm[l] * (
                            -real(Yb[i, l]) * sin(θ[i] - θ[l])
                            +
                            imag(Yb[i, l]) * cos(θ[i] - θ[l])
                        ) for l in bus_neighbors if l != i
                    )
                    Qi_Viθi = sum(
                        Vm[l] * (
                            real(Yb[i, l]) * cos(θ[i] - θ[l])
                            +
                            imag(Yb[i, l]) * sin(θ[i] - θ[l])
                        ) for l in bus_neighbors if l != i
                    )
                    Pi_ViVi = 2 * real(Yb[i, i])
                    Qi_ViVi = -2 * imag(Yb[i, i])

                    ViVis = Pi_ViVi * F_value[2 * i - 1] + Qi_ViVi * F_value[2 * i]
                    Viθis = Pi_Viθi * F_value[2 * i - 1] + Qi_Viθi * F_value[2 * i]

                    Hv[2 * i - 1, 2 * i - 1] += ViVis
                    Hv[2 * i, 2 * i - 1] += Viθis
                    Hv[2 * i - 1, 2 * i] += Viθis
                end
            else
                bt_k = data.bus_type[k, time_step]
                Gik, Bik = real(Yb[i, k]), imag(Yb[i, k])
                has_θk = (bt_k == PSY.ACBusTypes.PQ) || (bt_k == PSY.ACBusTypes.PV)
                has_θi = (bt_i == PSY.ACBusTypes.PQ) || (bt_i == PSY.ACBusTypes.PV)
                # the partials where all 3 indices are different vanish
                # naively count 8 with 2 distinct indices: {∂Vₖ, ∂θₖ} x {∂Vₖ, ∂θₖ, ∂Vᵢ, ∂θᵢ}
                # but can reduce to 6: ∂²/∂Vₖ∂θₖ = ∂²/∂θₖ∂Vₖ, and ∂²Δ{Pᵢ, Qᵢ}/∂²Vₖ is 0.
                # start with the 4 involving ∂θₖ, then do remaining the 2 involving ∂Vₖ
                if has_θk
                    # ∂²Δ{Pᵢ, Qᵢ}/∂²θₖ
                    Pi_θkθk =
                        Vm[i] * Vm[k] * ( # = Vm[k] * Qi_θkVk
                            -Gik * cos(θ[i] - θ[k])
                            -
                            Bik * sin(θ[i] - θ[k])
                        )
                    Qi_θkθk =
                        Vm[i] * Vm[k] * ( # = -Vm[k] * Pi_θkVk
                            -Gik * sin(θ[i] - θ[k])
                            +
                            Bik * cos(θ[i] - θ[k])
                        )
                    θkθks = Pi_θkθk * F_value[2 * i - 1] + Qi_θkθk * F_value[2 * i]
                    Hv[2 * k, 2 * k] += θkθks
                end
                if bt_k == PSY.ACBusTypes.PQ
                    # ∂²Δ{Pᵢ, Qᵢ}/∂θₖ∂Vₖ
                    Pi_θkVk = Vm[i] * (
                        Gik * sin(θ[i] - θ[k])
                        -
                        Bik * cos(θ[i] - θ[k])
                    )
                    Qi_θkVk = Vm[i] * (
                        -Gik * cos(θ[i] - θ[k])
                        -
                        Bik * sin(θ[i] - θ[k])
                    )
                    θkVks = Pi_θkVk * F_value[2 * i - 1] + Qi_θkVk * F_value[2 * i]
                    Hv[2 * k - 1, 2 * k] += θkVks
                    Hv[2 * k, 2 * k - 1] += θkVks
                end
                if has_θk && has_θi
                    # ∂²Δ{Pᵢ, Qᵢ}/∂θₖ∂θᵢ
                    Pi_θkθi =
                        Vm[i] * Vm[k] * (
                            Gik * cos(θ[i] - θ[k]) +
                            Bik * sin(θ[i] - θ[k])
                        )
                    Qi_θkθi =
                        Vm[i] * Vm[k] * (
                            Gik * sin(θ[i] - θ[k])
                            -
                            Bik * cos(θ[i] - θ[k])
                        )
                    θiθks = Pi_θkθi * F_value[2 * i - 1] + Qi_θkθi * F_value[2 * i]
                    Hv[2 * i, 2 * k] += θiθks
                    Hv[2 * k, 2 * i] += θiθks
                end
                if bt_i == PSY.ACBusTypes.PQ && has_θk
                    # ∂²Δ{Pᵢ, Qᵢ}/∂θₖ∂Vᵢ
                    Pi_θkVi = Vm[k] * ( # = Vm[k] * Qi_VkVi 
                        Gik * sin(θ[i] - θ[k])
                        -
                        Bik * cos(θ[i] - θ[k])
                    )
                    Qi_θkVi = Vm[k] * ( # = -Vm[k] * Pi_VkVi 
                        -Gik * cos(θ[i] - θ[k])
                        -
                        Bik * sin(θ[i] - θ[k])
                    )
                    Viθks = Pi_θkVi * F_value[2 * i - 1] + Qi_θkVi * F_value[2 * i]
                    Hv[2 * i - 1, 2 * k] += Viθks
                    Hv[2 * k, 2 * i - 1] += Viθks
                end
                if bt_k == PSY.ACBusTypes.PQ && has_θi
                    # ∂²Δ{Pᵢ, Qᵢ}/∂Vₖ∂θᵢ
                    Pi_Vkθi = Vm[i] * ( # = -Vm[i] * Qi_VkVi 
                        -Gik * sin(θ[i] - θ[k])
                        +
                        Bik * cos(θ[i] - θ[k])
                    )
                    Qi_Vkθi = Vm[i] * ( # = Vm[i] * Pi_VkVi 
                        Gik * cos(θ[i] - θ[k])
                        +
                        Bik * sin(θ[i] - θ[k])
                    )
                    θiVks = Pi_Vkθi * F_value[2 * i - 1] + Qi_Vkθi * F_value[2 * i]
                    Hv[2 * i, 2 * k - 1] += θiVks
                    Hv[2 * k - 1, 2 * i] += θiVks
                end
                if bt_k == PSY.ACBusTypes.PQ && bt_i == PSY.ACBusTypes.PQ
                    # ∂²Δ{Pᵢ, Qᵢ}/∂Vₖ∂Vᵢ
                    Pi_VkVi = Gik * cos(θ[i] - θ[k]) + Bik * sin(θ[i] - θ[k])
                    Qi_VkVi = Gik * sin(θ[i] - θ[k]) - Bik * cos(θ[i] - θ[k])
                    ViVks = Pi_VkVi * F_value[2 * i - 1] + Qi_VkVi * F_value[2 * i]
                    Hv[2 * i - 1, 2 * k - 1] += ViVks
                    Hv[2 * k - 1, 2 * i - 1] += ViVks
                end
            end
        end
    end
    return
end
