const ACPowerFlowData = PowerFlowData{
    PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
    Nothing,
}

struct ACPowerFlowResidual
    data::ACPowerFlowData
    Rf!::Function
    Rv::Vector{Float64}
    P_net::Vector{Float64}
    Q_net::Vector{Float64}
end

function ACPowerFlowResidual(data::ACPowerFlowData, time_step::Int64)
    n_buses = first(size(data.bus_type))
    P_net = zeros(n_buses)
    Q_net = zeros(n_buses)
    for ix in 1:n_buses
        P_net[ix] =
            data.bus_activepower_injection[ix, time_step] -
            data.bus_activepower_withdrawals[ix, time_step]
        Q_net[ix] =
            data.bus_reactivepower_injection[ix, time_step] -
            data.bus_reactivepower_withdrawals[ix, time_step]
    end

    return ACPowerFlowResidual(
        data,
        _update_residual_values!,
        Vector{Float64}(undef, 2 * n_buses),
        P_net,
        Q_net,
    )
end

function (Residual::ACPowerFlowResidual)(
    Rv::Vector{Float64},
    x::Vector{Float64},
    time_step::Int64,
)
    Residual.Rf!(Residual.Rv, x, Residual.P_net, Residual.Q_net, Residual.data, time_step)
    copyto!(Rv, Residual.Rv)
    return
end

function (Residual::ACPowerFlowResidual)(x::Vector{Float64}, time_step::Int64)
    Residual.Rf!(Residual.Rv, x, Residual.P_net, Residual.Q_net, Residual.data, time_step)
    return
end

# dispatching on Val for performance reasons.
function _set_state_vars_at_bus(
    ix::Int,
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.REF})
    # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
    P_net[ix] = StateVector[2 * ix - 1]
    Q_net[ix] = StateVector[2 * ix]
    data.bus_activepower_injection[ix, time_step] =
        StateVector[2 * ix - 1] + data.bus_activepower_withdrawals[ix, time_step]
    data.bus_reactivepower_injection[ix, time_step] =
        StateVector[2 * ix] + data.bus_reactivepower_withdrawals[ix, time_step]
end

function _set_state_vars_at_bus(
    ix::Int,
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.PV})
    # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
    Q_net[ix] = StateVector[2 * ix - 1]
    data.bus_reactivepower_injection[ix, time_step] =
        StateVector[2 * ix - 1] + data.bus_reactivepower_withdrawals[ix, time_step]
    data.bus_angles[ix, time_step] = StateVector[2 * ix]
end

function _set_state_vars_at_bus(
    ix::Int,
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.PQ})
    # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
    data.bus_magnitude[ix, time_step] = StateVector[2 * ix - 1]
    data.bus_angles[ix, time_step] = StateVector[2 * ix]
end

function _update_residual_values!(
    F::Vector{Float64},
    x::Vector{Float64},
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
)
    # update P_net, Q_net, data.bus_angles, data.bus_magnitude based on X
    Yb = data.power_network_matrix.data
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
        _set_state_vars_at_bus(ix, P_net, Q_net, x, data, time_step, Val(bt))
    end

    # compute active, reactive power balances using the just updated values.
    Vm = view(data.bus_magnitude, :, time_step)
    θ = view(data.bus_angles, :, time_step)
    # F is active and reactive power balance equations at all buses
    for bus_from in eachindex(P_net)
        S_re = 0.0
        S_im = 0.0
        for bus_to in data.neighbors[bus_from]
            gb = real(Yb[bus_from, bus_to])
            bb = imag(Yb[bus_from, bus_to])
            if bus_from == bus_to
                S_re += Vm[bus_from] * Vm[bus_to] * gb
                S_im += -Vm[bus_from] * Vm[bus_to] * bb
            else
                S_re +=
                    Vm[bus_from] *
                    Vm[bus_to] *
                    (gb * cos(θ[bus_from] - θ[bus_to]) + bb * sin(θ[bus_from] - θ[bus_to]))
                S_im +=
                    Vm[bus_from] *
                    Vm[bus_to] *
                    (gb * sin(θ[bus_from] - θ[bus_to]) - bb * cos(θ[bus_from] - θ[bus_to]))
            end
        end
        F[2 * bus_from - 1] = S_re - P_net[bus_from]
        F[2 * bus_from] = S_im - Q_net[bus_from]
    end
    return
end
