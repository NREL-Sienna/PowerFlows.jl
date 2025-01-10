
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
    return res.f_converged, res.zero
end