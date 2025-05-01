"""Partitions the state vector's variables based on what physical quantity each represents. 
Returns a `NamedTuple`, with the 4 keys `Va`, `Vm`, `P`, and `Q`. The 4 values are vectors 
of length equal to the number of buses, with `NaN`s in the positions where that physical 
quantity is not part of the state vector for that bus. (Currently not intended for use in 
spots where performance is critical.)"""
function partition_state(x::Vector{Float64},
    bus_types::AbstractVector{PSY.ACBusTypes},
)
    # usually, bus_types will be data.bus_type[:, time_step]
    nbuses = div(size(x, 1), 2)
    (Vms, Vas, Ps, Qs) = (fill(NaN, nbuses) for _ in 1:4)
    for i in 1:nbuses
        if bus_types[i] == PSY.ACBusTypes.REF
            Ps[i] = x[2 * i - 1]
            Qs[i] = x[2 * i]
        elseif bus_types[i] == PSY.ACBusTypes.PQ
            Vms[i] = x[2 * i - 1]
            Vas[i] = x[2 * i]
        elseif bus_types[i] == PSY.ACBusTypes.PV
            Qs[i] = x[2 * i - 1]
            Vas[i] = x[2 * i]
        end
    end
    return (Va = Vas, Vm = Vms, P = Ps, Q = Qs)
end

"""Update state vector based on values of fields of data."""
function update_state!(x::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
)
    # I really only need access to 3 classes of fields of data:
    # bus types, bus power contributions [inj/widthdrawal], and bus voltages
    @assert length(x) == 2 * length(data.bus_type[:, 1])
    state_variable_count = 1
    for (ix, b) in enumerate(data.bus_type[:, time_step])
        if b == PSY.ACBusTypes.REF
            x[state_variable_count] =
                data.bus_activepower_injection[ix, time_step] -
                data.bus_activepower_withdrawals[ix, time_step]
            x[state_variable_count + 1] =
                data.bus_reactivepower_injection[ix, time_step] -
                data.bus_reactivepower_withdrawals[ix, time_step]
            state_variable_count += 2
        elseif b == PSY.ACBusTypes.PV
            x[state_variable_count] =
                data.bus_reactivepower_injection[ix, time_step] -
                data.bus_reactivepower_withdrawals[ix, time_step]
            x[state_variable_count + 1] = data.bus_angles[ix, time_step]
            state_variable_count += 2
        elseif b == PSY.ACBusTypes.PQ
            x[state_variable_count] = data.bus_magnitude[ix, time_step]
            x[state_variable_count + 1] = data.bus_angles[ix, time_step]
            state_variable_count += 2
        else
            throw(ArgumentError("$b not recognized as a bustype"))
        end
    end
    @assert state_variable_count - 1 == length(data.bus_type[:, 1]) * 2
end

"""Update the fields of data based on the values of the state vector."""
function update_data!(data::ACPowerFlowData,
    x::Vector{Float64},
    time_step::Int64,
)
    # same as above. I really only need access to 3 classes of fields of data:
    # bus types, bus power contributions [inj/widthdrawal], and bus voltage [angle/mag]
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
        _set_state_vars_at_bus(ix, x, data, time_step, Val(bt))
    end
end

# dispatching on Val for performance reasons.
function _set_state_vars_at_bus(
    ix::Int,
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.REF})
    # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
    data.bus_activepower_injection[ix, time_step] =
        StateVector[2 * ix - 1] + data.bus_activepower_withdrawals[ix, time_step]
    data.bus_reactivepower_injection[ix, time_step] =
        StateVector[2 * ix] + data.bus_reactivepower_withdrawals[ix, time_step]
end

function _set_state_vars_at_bus(
    ix::Int,
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.PV})
    # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
    data.bus_reactivepower_injection[ix, time_step] =
        StateVector[2 * ix - 1] + data.bus_reactivepower_withdrawals[ix, time_step]
    data.bus_angles[ix, time_step] = StateVector[2 * ix]
end

function _set_state_vars_at_bus(
    ix::Int,
    x::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.PQ})
    # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
    data.bus_magnitude[ix, time_step] = x[2 * ix - 1]
    data.bus_angles[ix, time_step] = x[2 * ix]
end

# could dispatch on val, and/or combine with update_data!.
"""Update `P_net` and `Q_net` based on the values of the state vector."""
function update_net_power!(P_net::Vector{Float64},
    Q_net::Vector{Float64},
    x::Vector{Float64},
    bus_types::AbstractVector{PSY.ACBusTypes},
)
    # usually, bus_types will be data.bus_type[:, time_step]
    for (ix, bt) in enumerate(bus_types)
        if bt == PSY.ACBusTypes.REF
            P_net[ix] = x[2 * ix - 1]
            Q_net[ix] = x[2 * ix]
        elseif bt == PSY.ACBusTypes.PV
            Q_net[ix] = x[2 * ix - 1]
        end
    end
    return
end
