const ACPowerFlowData = PowerFlowData{
    PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
    Nothing,
    String,
}

struct PolarPowerFlow{F, D}
    Pf::F
    data::D
    residual::Vector{Float64}
    P_net::Vector{Float64}
    Q_net::Vector{Float64}
end

function PolarPowerFlow(data::ACPowerFlowData)
    n_buses = length(data.bus_type)
    P_net = zeros(n_buses)
    Q_net = zeros(n_buses)
    for ix in 1:n_buses
        P_net[ix] =
            data.bus_activepower_injection[ix] - data.bus_activepower_withdrawals[ix]
        Q_net[ix] =
            data.bus_reactivepower_injection[ix] - data.bus_reactivepower_withdrawals[ix]
    end
    return PolarPowerFlow(
        (res::Vector{Float64}, X::Vector{Float64}) -> polar_pf!(res, X, P_net, Q_net, data),
        data,
        zeros(2 * n_buses),
        P_net,
        Q_net,
    )
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
    data::ACPowerFlowData,
)
    Yb = data.power_network_matrix
    n_buses = length(data.bus_type)
    for (ix, b) in enumerate(data.bus_type)
        if b == PSY.ACBusTypes.REF
            # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
            P_net[ix] = X[2 * ix - 1] - data.bus_activepower_withdrawals[ix]
            Q_net[ix] = X[2 * ix] - data.bus_reactivepower_withdrawals[ix]
        elseif b == PSY.ACBusTypes.PV
            # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
            Q_net[ix] = X[2 * ix - 1] - data.bus_reactivepower_withdrawals[ix]
            data.bus_angles[ix] = X[2 * ix]
        elseif b == PSY.ACBusTypes.PQ
            # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
            data.bus_magnitude[ix] = X[2 * ix - 1]
            data.bus_angles[ix] = X[2 * ix]
        end
    end

    Vm = data.bus_magnitude
    θ = data.bus_angles
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
