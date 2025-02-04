const ACPowerFlowData = PowerFlowData{
    PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
    Nothing,
}

struct PolarPowerFlow{F, D}
    Pf::F
    data::D
    residual::Vector{Float64}
    P_net::Vector{Float64}
    Q_net::Vector{Float64}
    x0::Vector{Float64}
end

function _calculate_x0(time_step::Int64,
    bus_types::Matrix{PSY.ACBusTypes},
    bus_angles::Matrix{Float64},
    bus_magnitude::Matrix{Float64},
    bus_activepower_injection::Matrix{Float64},
    bus_reactivepower_injection::Matrix{Float64},
    bus_activepower_withdrawals::Matrix{Float64},
    bus_reactivepower_withdrawals::Matrix{Float64})
    n_buses = length(bus_types[:, 1])
    x0 = Vector{Float64}(undef, 2 * n_buses)
    state_variable_count = 1
    for (ix, b) in enumerate(bus_types[:, time_step])
        if b == PSY.ACBusTypes.REF
            x0[state_variable_count] =
                bus_activepower_injection[ix, time_step] -
                bus_activepower_withdrawals[ix, time_step]
            x0[state_variable_count + 1] =
                bus_reactivepower_injection[ix, time_step] -
                bus_reactivepower_withdrawals[ix, time_step]
            state_variable_count += 2
        elseif b == PSY.ACBusTypes.PV
            x0[state_variable_count] =
                bus_reactivepower_injection[ix, time_step] -
                bus_reactivepower_withdrawals[ix, time_step]
            x0[state_variable_count + 1] = bus_angles[ix, time_step]
            state_variable_count += 2
        elseif b == PSY.ACBusTypes.PQ
            x0[state_variable_count] = bus_magnitude[ix, time_step]
            x0[state_variable_count + 1] = bus_angles[ix, time_step]
            state_variable_count += 2
        else
            throw(ArgumentError("$b not recognized as a bustype"))
        end
    end
    @assert state_variable_count - 1 == n_buses * 2
    return x0
end

function PolarPowerFlow(data::ACPowerFlowData; time_step::Int64 = 1)
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
    x0 = _calculate_x0(time_step,
        data.bus_type,
        data.bus_angles,
        data.bus_magnitude,
        data.bus_activepower_injection,
        data.bus_reactivepower_injection,
        data.bus_activepower_withdrawals,
        data.bus_reactivepower_withdrawals)
    pf_function =
        (res::Vector{Float64}, X::Vector{Float64}) ->
            polar_pf!(res, X, P_net, Q_net, data; time_step = time_step)
    res = similar(x0)
    pf_function(res, x0)

    if sum(res) > 10 * (n_buses * 2)
        _, ix = findmax(res)
        bus_no = data.bus_lookup[ix]
        @warn "Initial guess provided results in a large initial residual. Largest residual at bus $bus_no"
    end

    return PolarPowerFlow(
        pf_function,
        data,
        res,
        P_net,
        Q_net,
        x0,
    )
end

function (pf::PolarPowerFlow)(res::Vector{Float64}, x::Vector{Float64})
    pf.Pf(pf.residual, x)
    copyto!(res, pf.residual)
    return
end

function (pf::PolarPowerFlow)(x::Vector{Float64})
    pf.Pf(pf.residual, x)
    return
end

function polar_pf!(
    F::Vector{Float64},
    X::Vector{Float64},
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    data::ACPowerFlowData;
    time_step::Int64 = 1,
)
    Yb = data.power_network_matrix.data
    n_buses = first(size(data.bus_type))
    for (ix, b) in enumerate(data.bus_type[:, time_step])
        if b == PSY.ACBusTypes.REF
            # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
            P_net[ix] = X[2 * ix - 1] - data.bus_activepower_withdrawals[ix, time_step]
            Q_net[ix] = X[2 * ix] - data.bus_reactivepower_withdrawals[ix, time_step]
        elseif b == PSY.ACBusTypes.PV
            # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
            Q_net[ix] = X[2 * ix - 1] - data.bus_reactivepower_withdrawals[ix, time_step]
            data.bus_angles[ix, time_step] = X[2 * ix]
        elseif b == PSY.ACBusTypes.PQ
            # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
            data.bus_magnitude[ix, time_step] = X[2 * ix - 1]
            data.bus_angles[ix, time_step] = X[2 * ix]
        end
    end

    Vm = data.bus_magnitude[:, time_step]
    θ = data.bus_angles[:, time_step]
    # F is active and reactive power balance equations at all buses
    for ix_f in 1:n_buses
        S_re = -P_net[ix_f]
        S_im = -Q_net[ix_f]
        for ix_t in data.neighbors[ix_f]
            gb = real(Yb[ix_f, ix_t])
            bb = imag(Yb[ix_f, ix_t])
            if ix_f == ix_t
                S_re += Vm[ix_f] * Vm[ix_t] * gb
                S_im += -Vm[ix_f] * Vm[ix_t] * bb
            else
                S_re +=
                    Vm[ix_f] *
                    Vm[ix_t] *
                    (gb * cos(θ[ix_f] - θ[ix_t]) + bb * sin(θ[ix_f] - θ[ix_t]))
                S_im +=
                    Vm[ix_f] *
                    Vm[ix_t] *
                    (gb * sin(θ[ix_f] - θ[ix_t]) - bb * cos(θ[ix_f] - θ[ix_t]))
            end
        end
        F[2 * ix_f - 1] = S_re
        F[2 * ix_f] = S_im
    end
    return
end
