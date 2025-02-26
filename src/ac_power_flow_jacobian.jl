struct PolarPowerFlowJacobian{F,
    T <: SparseArrays.SparseMatrixCSC{Float64, Int32},
}
    Jf::F
    Jv::T
end

"""Update the Jacobian. Point of confusion: Jf is jsp!(), which
updates Jv based on data, not based on x. x isn't actually used."""
function (Jacobianf::PolarPowerFlowJacobian)(
    x::Vector{Float64},
)
    Jacobianf.Jf(Jacobianf.Jv, x)
    return
end

"""Update the Jacobian, and write new Jacobian to J. Unused currently."""
function (Jacobianf::PolarPowerFlowJacobian)(
    J::SparseArrays.SparseMatrixCSC{Float64, Int32},
    x::Vector{Float64},
)
    Jacobianf.Jf(Jacobianf.Jv, x)
    copyto!(J, Jacobianf.Jv)
    return
end

function PolarPowerFlowJacobian(
    data::ACPowerFlowData,
    x0::Vector{Float64},
    time_step::Int64,
)
    j0 = _create_jacobian_matrix_structure(data, time_step)
    # yay closures: changes to data will affect F0 and vice versa.
    F0 =
        (J::SparseArrays.SparseMatrixCSC{Float64, Int32}, x::Vector{Float64}) ->
            jsp!(J, x, data, time_step)
    F0(j0, x0)
    return PolarPowerFlowJacobian(F0, j0)
end

function _create_jacobian_matrix_structure(data::ACPowerFlowData, time_step::Int64)
    #Create Jacobian structure
    J0_I = Int32[]
    J0_J = Int32[]
    J0_V = Float64[]

    num_buses = first(size(data.bus_type))

    for ix_f in 1:num_buses
        F_ix_f_r = 2 * ix_f - 1
        F_ix_f_i = 2 * ix_f
        for ix_t in data.neighbors[ix_f]
            X_ix_t_fst = 2 * ix_t - 1
            X_ix_t_snd = 2 * ix_t
            nb = data.bus_type[ix_t, time_step]
            #Set to 0.0 only on connected buses
            if nb == PSY.ACBusTypes.REF
                if ix_f == ix_t
                    #Active PF w/r Local Active Power
                    push!(J0_I, F_ix_f_r)
                    push!(J0_J, X_ix_t_fst)
                    push!(J0_V, 0.0)
                    #Rective PF w/r Local Reactive Power
                    push!(J0_I, F_ix_f_i)
                    push!(J0_J, X_ix_t_snd)
                    push!(J0_V, 0.0)
                end
            elseif nb == PSY.ACBusTypes.PV
                #Active PF w/r Angle
                push!(J0_I, F_ix_f_r)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
                #Reactive PF w/r Angle
                push!(J0_I, F_ix_f_i)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
                if ix_f == ix_t
                    #Reactive PF w/r Local Reactive Power
                    push!(J0_I, F_ix_f_i)
                    push!(J0_J, X_ix_t_fst)
                    push!(J0_V, 0.0)
                end
            elseif nb == PSY.ACBusTypes.PQ
                #Active PF w/r VoltageMag
                push!(J0_I, F_ix_f_r)
                push!(J0_J, X_ix_t_fst)
                push!(J0_V, 0.0)
                #Reactive PF w/r VoltageMag
                push!(J0_I, F_ix_f_i)
                push!(J0_J, X_ix_t_fst)
                push!(J0_V, 0.0)
                #Active PF w/r Angle
                push!(J0_I, F_ix_f_r)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
                #Reactive PF w/r Angle
                push!(J0_I, F_ix_f_i)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
            end
        end
    end
    J0 = SparseArrays.sparse(J0_I, J0_J, J0_V)
    return J0
end

function _set_entries_for_neighbor(J::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Yb::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    Vm::Vector{Float64},
    θ::Vector{Float64},
    ix_f::Int,
    ix_t::Int,
    F_ix_f_r::Int,
    F_ix_f_i::Int,
    X_ix_t_fst::Int,
    X_ix_t_snd::Int,
    ix_f_neighbors::Set{Int},
    ::Val{PSY.ACBusTypes.REF})
    # State variables are Active and Reactive Power Generated
    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - x[2*i-1]
    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
    # x does not appears in p_flow and q_flow
    if ix_f == ix_t
        J[F_ix_f_r, X_ix_t_fst] = -1.0
        J[F_ix_f_i, X_ix_t_snd] = -1.0
    end
    return
end

function _set_entries_for_neighbor(J::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Yb::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    Vm::Vector{Float64},
    θ::Vector{Float64},
    ix_f::Int,
    ix_t::Int,
    F_ix_f_r::Int,
    F_ix_f_i::Int,
    X_ix_t_fst::Int,
    X_ix_t_snd::Int,
    ix_f_neighbors::Set{Int},
    ::Val{PSY.ACBusTypes.PV})
    # State variables are Reactive Power Generated and Voltage Angle
    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - p_gen[i]
    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
    # x[2*i] (associated with q_gen) does not appear in q_flow
    # x[2*i] (associated with q_gen) does not appear in the active power balance
    if ix_f == ix_t
        #Jac: Reactive PF against local active power
        J[F_ix_f_i, X_ix_t_fst] = -1.0
        #Jac: Active PF against same Angle: θ[ix_f] =  θ[ix_t]
        J[F_ix_f_r, X_ix_t_snd] =
            Vm[ix_f] * sum(
                Vm[k] * (
                    real(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k]) +
                    imag(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k])
                ) for k in ix_f_neighbors if k != ix_f
            )
        #Jac: Reactive PF against same Angle: θ[ix_f] = θ[ix_t]
        J[F_ix_f_i, X_ix_t_snd] =
            Vm[ix_f] * sum(
                Vm[k] * (
                    real(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k]) -
                    imag(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k])
                ) for k in ix_f_neighbors if k != ix_f
            )
    else
        g_ij = real(Yb[ix_f, ix_t])
        b_ij = imag(Yb[ix_f, ix_t])
        #Jac: Active PF against other angles θ[ix_t]
        J[F_ix_f_r, X_ix_t_snd] =
            Vm[ix_f] *
            Vm[ix_t] *
            (g_ij * sin(θ[ix_f] - θ[ix_t]) + b_ij * -cos(θ[ix_f] - θ[ix_t]))
        #Jac: Reactive PF against other angles θ[ix_t]
        J[F_ix_f_i, X_ix_t_snd] =
            Vm[ix_f] *
            Vm[ix_t] *
            (g_ij * -cos(θ[ix_f] - θ[ix_t]) - b_ij * sin(θ[ix_f] - θ[ix_t]))
    end
    return
end

function _set_entries_for_neighbor(J::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Yb::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    Vm::Vector{Float64},
    θ::Vector{Float64},
    ix_f::Int,
    ix_t::Int,
    F_ix_f_r::Int,
    F_ix_f_i::Int,
    X_ix_t_fst::Int,
    X_ix_t_snd::Int,
    ix_f_neighbors::Set{Int},
    ::Val{PSY.ACBusTypes.PQ})
    # State variables are Voltage Magnitude and Voltage Angle
    # Everything appears in everything
    if ix_f == ix_t
        #Jac: Active PF against same voltage magnitude Vm[ix_f]
        J[F_ix_f_r, X_ix_t_fst] =
            2 * real(Yb[ix_f, ix_t]) * Vm[ix_f] + sum(
                Vm[k] * (
                    real(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k]) +
                    imag(Yb[ix_f, k]) * sin(θ[ix_f] - θ[k])
                ) for k in ix_f_neighbors if k != ix_f
            )
        #Jac: Active PF against same angle θ[ix_f]
        J[F_ix_f_r, X_ix_t_snd] =
            Vm[ix_f] * sum(
                Vm[k] * (
                    real(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k]) +
                    imag(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k])
                ) for k in ix_f_neighbors if k != ix_f
            )

        #Jac: Reactive PF against same voltage magnitude Vm[ix_f]
        J[F_ix_f_i, X_ix_t_fst] =
            -2 * imag(Yb[ix_f, ix_t]) * Vm[ix_f] + sum(
                Vm[k] * (
                    real(Yb[ix_f, k]) * sin(θ[ix_f] - θ[k]) -
                    imag(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k])
                ) for k in ix_f_neighbors if k != ix_f
            )
        #Jac: Reactive PF against same angle θ[ix_f]
        J[F_ix_f_i, X_ix_t_snd] =
            Vm[ix_f] * sum(
                Vm[k] * (
                    real(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k]) -
                    imag(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k])
                ) for k in ix_f_neighbors if k != ix_f
            )
    else
        g_ij = real(Yb[ix_f, ix_t])
        b_ij = imag(Yb[ix_f, ix_t])
        #Jac: Active PF w/r to different voltage magnitude Vm[ix_t]
        J[F_ix_f_r, X_ix_t_fst] =
            Vm[ix_f] *
            (g_ij * cos(θ[ix_f] - θ[ix_t]) + b_ij * sin(θ[ix_f] - θ[ix_t]))
        #Jac: Active PF w/r to different angle θ[ix_t]
        J[F_ix_f_r, X_ix_t_snd] =
            Vm[ix_f] *
            Vm[ix_t] *
            (g_ij * sin(θ[ix_f] - θ[ix_t]) + b_ij * -cos(θ[ix_f] - θ[ix_t]))

        #Jac: Reactive PF w/r to different voltage magnitude Vm[ix_t]
        J[F_ix_f_i, X_ix_t_fst] =
            Vm[ix_f] *
            (g_ij * sin(θ[ix_f] - θ[ix_t]) - b_ij * cos(θ[ix_f] - θ[ix_t]))
        #Jac: Reactive PF w/r to different angle θ[ix_t]
        J[F_ix_f_i, X_ix_t_snd] =
            Vm[ix_f] *
            Vm[ix_t] *
            (g_ij * -cos(θ[ix_f] - θ[ix_t]) - b_ij * sin(θ[ix_f] - θ[ix_t]))
    end
    return
end

"""Used to update Jv based on the bus voltages, angles, etc. in data."""
function jsp!(
    J::SparseArrays.SparseMatrixCSC{Float64, Int32},
    ::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
)
    Yb = data.power_network_matrix.data
    Vm = data.bus_magnitude[:, time_step]
    θ = data.bus_angles[:, time_step]
    num_buses = first(size(data.bus_type))
    for ix_f in 1:num_buses
        F_ix_f_r = 2 * ix_f - 1
        F_ix_f_i = 2 * ix_f

        for ix_t in data.neighbors[ix_f]
            X_ix_t_fst = 2 * ix_t - 1
            X_ix_t_snd = 2 * ix_t
            nb = data.bus_type[ix_t, time_step]
            # could move the definitions of 
            # F_ix_f_{i/r}, X_ix_t_{fst/snd} inside _set_entries.
            _set_entries_for_neighbor(J, Yb, Vm, θ,
                ix_f, ix_t,
                F_ix_f_r, F_ix_f_i,
                X_ix_t_fst, X_ix_t_snd,
                data.neighbors[ix_f], Val(nb))
        end
    end
    return
end

function calculate_loss_factors(
    data::ACPowerFlowData,
    Jv::SparseMatrixCSC{Float64, Int32},
    time_step::Int,
)
    num_buses = first(size(data.bus_type))
    ref, pv, pq = bus_type_idx(data, time_step)
    pvpq = vcat(pv, pq)
    pvpq_coords = [
        x for pair in zip(
            [2 * x - 1 for x in 1:num_buses if x in pvpq],
            [2 * x for x in 1:num_buses if x in pvpq],
        ) for x in pair
    ]
    data.loss_factors[ref, time_step] .= 0.0
    penalty_factors!(
        Jv[pvpq_coords, pvpq_coords],
        vec(collect(Jv[2 .* ref .- 1, pvpq_coords])),
        view(
            data.loss_factors,
            [x for x in 1:num_buses if x in pvpq],
            time_step,
        ),
        [2 * x - 1 for x in 1:length(pvpq)],
    )
end
