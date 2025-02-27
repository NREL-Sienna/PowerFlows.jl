const ACPowerFlowData = PowerFlowData{
    PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
    Nothing,
}

struct ACPowerFlowResidual{F}
    Pf::F
    data::ACPowerFlowData
    residual::Vector{Float64}
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
    pf_function =
        (res::Vector{Float64}, X::Vector{Float64}) ->
            polar_pf!(res, X, P_net, Q_net, data, time_step)

    return ACPowerFlowResidual(
        pf_function,
        data,
        Vector{Float64}(undef, 2*n_buses),
        P_net,
        Q_net,
    )
end

function (pf::ACPowerFlowResidual)(res::Vector{Float64}, x::Vector{Float64})
    pf.Pf(pf.residual, x)
    copyto!(res, pf.residual)
    return
end

function (pf::ACPowerFlowResidual)(x::Vector{Float64})
    pf.Pf(pf.residual, x)
    return
end

function polar_pf!(
    F::Vector{Float64},
    X::Vector{Float64},
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
)
    # update P_net, Q_net, data.bus_angles, data.bus_magnitude based on X
    Yb = data.power_network_matrix.data
    n_buses = first(size(data.bus_type))
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
        if bt == PSY.ACBusTypes.REF
            # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
            P_net[ix] = X[2 * ix - 1] - data.bus_activepower_withdrawals[ix, time_step]
            Q_net[ix] = X[2 * ix] - data.bus_reactivepower_withdrawals[ix, time_step]
        elseif bt == PSY.ACBusTypes.PV
            # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
            Q_net[ix] = X[2 * ix - 1] - data.bus_reactivepower_withdrawals[ix, time_step]
            data.bus_angles[ix, time_step] = X[2 * ix]
        elseif bt == PSY.ACBusTypes.PQ
            # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
            data.bus_magnitude[ix, time_step] = X[2 * ix - 1]
            data.bus_angles[ix, time_step] = X[2 * ix]
        end
    end

    # compute active, reactive power balances using the just updated values.
    Vm = data.bus_magnitude[:, time_step]
    θ = data.bus_angles[:, time_step]
    # F is active and reactive power balance equations at all buses
    for ix_f in 1:n_buses
        S_re = 0.0
        S_im = 0.0
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
        F[2 * ix_f - 1] = S_re - P_net[ix_f]
        F[2 * ix_f] = S_im - Q_net[ix_f]
    end
    return
end
