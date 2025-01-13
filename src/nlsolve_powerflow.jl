
function _newton_powerflow(
    pf::ACPowerFlow{NLSolveACPowerFlow},
    data::ACPowerFlowData;
    time_step::Int64 = 1,  # not implemented for NLSolve and not used
    nlsolve_kwargs...,
)
    if time_step != 1
        throw(
            ArgumentError(
                "Multiperiod power flow not implemented for NLSolve AC power flow",
            ),
        )
    end

    pf = PolarPowerFlow(data)
    J = PowerFlows.PolarPowerFlowJacobian(data, pf.x0)

    df = NLsolve.OnceDifferentiable(pf, J, pf.x0, pf.residual, J.Jv)
    res = NLsolve.nlsolve(df, pf.x0; nlsolve_kwargs...)
    if !res.f_converged
        @error(
            "The powerflow solver NLSolve did not converge (returned convergence = $(res.f_converged))"
        )
    end
    V = _calc_V(data, res.zero; time_step = time_step)
    Sbus_result = V .* conj(data.power_network_matrix.data * V)
    return (res.f_converged, V, Sbus_result)
end

function _calc_V(
    data::ACPowerFlowData,
    x::Vector{Float64};
    time_step::Int64 = 1,
)
    n_buses = length(x) ÷ 2  # Since x has 2 elements per bus (real and imaginary)
    V = zeros(Complex{Float64}, n_buses)
    Vm_data = data.bus_magnitude[:, time_step]
    Va_data = data.bus_angles[:, time_step]

    # Extract values for Vm and Va from x
    for (ix, b) in enumerate(data.bus_type)
        if b == PSY.ACBusTypes.REF
            # For REF bus, we have active and reactive power
            Vm = Vm_data[ix]
            Va = Va_data[ix]
            V[ix] = Vm * exp(im * Va)
        elseif b == PSY.ACBusTypes.PV
            # For PV bus, we have reactive power and voltage angle
            Vm = Vm_data[ix]
            Va = x[2 * ix]
            V[ix] = Vm * exp(im * Va)  # Rebuild voltage from magnitude and angle
        elseif b == PSY.ACBusTypes.PQ
            # For PQ bus, we have voltage magnitude and voltage angle
            Vm = x[2 * ix - 1]
            Va = x[2 * ix]
            V[ix] = Vm * exp(im * Va)  # Rebuild voltage from magnitude and angle
        end
    end

    return V
end
