struct HomotopyHessian
    data::ACPowerFlowData
    pfResidual::ACPowerFlowResidual
    J::ACPowerFlowJacobian
    Hv::SparseMatrixCSC{Float64, Int32}
    t_k_ref::Base.RefValue{Float64}
end

function (hess::HomotopyHessian)(Hv::SparseMatrixCSC{Float64, Int32},
    x::Vector{Float64};
    internal = false)
    println("hessian at $x")
    time_step = 1
    hess.pfResidual(x, time_step)
    hess.J(time_step)
    _update_hessian_matrix_values(Hv, hess.pfResidual.Rv, hess.data, time_step)
    Jv = hess.J.Jv
    Hv .+= Jv' * Jv
    SparseArrays.nonzeros(hess.Hv) .*= hess.t_k_ref[]
    for (bus_ix, bus_type) in enumerate(get_bus_type(hess.data)[:, time_step])
        if bus_type == PSY.ACBusTypes.PQ
            Hv[2 * bus_ix - 1, 2 * bus_ix - 1] += 1.0 - hess.t_k_ref[]
        end
    end
    if !internal
        copyto!(hess.Hv, Hv)
    end
    return nothing
end

function (hess::HomotopyHessian)(x::Vector{Float64})
    hess(hess.Hv, x; internal = true)
    return nothing
end

function HomotopyHessian(data::ACPowerFlowData, time_step::Int)
    pfResidual = ACPowerFlowResidual(data, time_step)
    Hv = _create_hessian_matrix_structure(data, time_step)
    J = ACPowerFlowJacobian(data, time_step)
    return HomotopyHessian(data, pfResidual, J, Hv, Ref(0.0))
end

function _create_hessian_matrix_structure(data::ACPowerFlowData, time_step::Int64)
    rows = Int32[]      # I
    columns = Int32[]   # J
    values = Float64[]  # V

    num_buses = first(size(data.bus_type))

    for bus_from in 1:num_buses
        # I'll follow the same strategy as the Jacobian, and take the worst case
        # where all buses are PQ.
        for bus_to in data.neighbors[bus_from]
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

"""Sets Hv equal to F_1(x) H_{F_1}(x) + ...+ F_{2n}(x) H_{F_{2n}}(x),
where F_k denotes the kth power balance equation and H_{F_k} its Hessian."""
function _update_hessian_matrix_values(
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
end
