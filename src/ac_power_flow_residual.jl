const ACPowerFlowData = PowerFlowData{
    PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
    Nothing,
}

"""
    struct ACPowerFlowResidual

A struct to keep track of the residuals in the Newton-Raphson AC power flow calculation.

# Fields
- `data::ACPowerFlowData`: The grid model data.
- `Rf!::Function`: A function that updates the residuals based on the latest values stored in the grid at the given iteration.
- `Rv::Vector{Float64}`: A vector of the values of the residuals.
- `P_net::Vector{Float64}`: A vector of net active power injections.
- `Q_net::Vector{Float64}`: A vector of net reactive power injections.
"""
struct ACPowerFlowResidual
    data::ACPowerFlowData
    Rf!::Function
    Rv::Vector{Float64}
    P_net::Vector{Float64}
    Q_net::Vector{Float64}
end

"""
    ACPowerFlowResidual(data::ACPowerFlowData, time_step::Int64)

Create an instance of `ACPowerFlowResidual` for a given time step.

# Arguments
- `data::ACPowerFlowData`: The power flow data representing the power system model.
- `time_step::Int64`: The time step for which the power flow calculation is executed.

# Returns
- `ACPowerFlowResidual`: An instance containing the residual values, net bus active power injections, 
and net bus reactive power injections.
"""
function ACPowerFlowResidual(data::ACPowerFlowData, time_step::Int64)
    n_buses = first(size(data.bus_type))
    P_net = Vector{Float64}(undef, n_buses)
    Q_net = Vector{Float64}(undef, n_buses)
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

"""
    (Residual::ACPowerFlowResidual)(Rv::Vector{Float64}, x::Vector{Float64}, time_step::Int64)

Evaluate the AC power flow residuals and store the result in `Rv` using the provided 
state vector `x` and the current time step `time_step`.
The residuals are updated inplace in the struct and additionally copied to the provided array. 
This function implements the functor approach for the `ACPowerFlowResidual` struct. 
This makes the struct callable.
Calling the `ACPowerFlowResidual` will also update the values of P, Q, V, Θ in the `data` struct.

# Arguments
- `Rv::Vector{Float64}`: The vector to store the calculated residuals.
- `x::Vector{Float64}`: The state vector.
- `time_step::Int64`: The current time step.
"""
function (Residual::ACPowerFlowResidual)(
    Rv::Vector{Float64},
    x::Vector{Float64},
    time_step::Int64,
)
    Residual.Rf!(Residual.Rv, x, Residual.P_net, Residual.Q_net, Residual.data, time_step)
    copyto!(Rv, Residual.Rv)
    return
end

"""
    (Residual::ACPowerFlowResidual)(x::Vector{Float64}, time_step::Int64)

Update the AC power flow residuals inplace and store the result in the attribute `Rv` of the struct.
The inputs are the values of state vector `x` and the current time step `time_step`.
This function implements the functor approach for the `ACPowerFlowResidual` struct. 
This makes the struct callable.
Calling the `ACPowerFlowResidual` will also update the values of P, Q, V, Θ in the `data` struct.

# Arguments
- `x::Vector{Float64}`: The state vector values.
- `time_step::Int64`: The current time step.
"""
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

"""
    _update_residual_values!(
        F::Vector{Float64},
        x::Vector{Float64},
        P_net::Vector{Float64},
        Q_net::Vector{Float64},
        data::ACPowerFlowData,
        time_step::Int64,
    )

Update the residual values for the Newton-Raphson AC power flow calculation. This function is used internally in the
`ACPowerFlowResidual` struct. This function also updates the values of P, Q, V, Θ in the `data` struct.

# Arguments
- `F::Vector{Float64}`: Vector of the values of the residuals.
- `x::Vector{Float64}`: State vector values.
- `P_net::Vector{Float64}`: Vector of net active power injections at each bus.
- `Q_net::Vector{Float64}`: Vector of net reactive power injections at each bus.
- `data::ACPowerFlowData`: Data structure representing the grid model for the AC power flow calculation.
- `time_step::Int64`: The current time step for which the residual values are being updated.
"""
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
