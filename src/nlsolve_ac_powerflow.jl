"""
Solves a the power flow into the system and writes the solution into the relevant structs.
Updates generators active and reactive power setpoints and branches active and reactive
power flows (calculated in the From - To direction) (see
[`flow_val`](@ref))

Supports passing NLsolve kwargs in the args. By default shows the solver trace.

Arguments available for `nlsolve`:

  - `get_connectivity::Bool`: Checks if the network is connected. Default true
  - `method` : See NLSolve.jl documentation for available solvers
  - `xtol`: norm difference in `x` between two successive iterates under which
    convergence is declared. Default: `0.0`.
  - `ftol`: infinite norm of residuals under which convergence is declared.
    Default: `1e-8`.
  - `iterations`: maximum number of iterations. Default: `1_000`.
  - `store_trace`: should a trace of the optimization algorithm's state be
    stored? Default: `false`.
  - `show_trace`: should a trace of the optimization algorithm's state be shown
    on `STDOUT`? Default: `false`.
  - `extended_trace`: should additifonal algorithm internals be added to the state
    trace? Default: `false`.

## Examples

```julia
solve_powerflow!(sys)
# Passing NLsolve arguments
solve_powerflow!(sys, method=:newton)
```
"""
function solve_ac_powerflow!(system::PSY.System; kwargs...)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    data = PowerFlowData(
        ACPowerFlow(),
        system;
        check_connectivity = get(kwargs, :check_connectivity, true),
    )
    res = _solve_powerflow(data; kwargs...)
    if res.f_converged
        write_powerflow_solution!(system, res.zero)
        @info("PowerFlow solve converged, the results have been stored in the system")
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return res.f_converged
    end
    @error("The powerflow solver returned convergence = $(res.f_converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return res.f_converged
end

"""
Similar to solve_powerflow!(sys) but does not update the system struct with results.
Returns the results in a dictionary of dataframes.

## Examples

```julia
res = solve_powerflow(sys)
# Passing NLsolve arguments
res = solve_powerflow(sys, method=:newton)
```
"""
function solve_powerflow(
    ::ACPowerFlow,
    system::PSY.System;
    kwargs...,
)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    data = PowerFlowData(
        ACPowerFlow(),
        system;
        check_connectivity = get(kwargs, :check_connectivity, true),
    )
    res = _solve_powerflow(data; kwargs...)
    if res.f_converged
        @info("PowerFlow solve converged, the results are exported in DataFrames")
        df_results = write_results(ACPowerFlow(), system, res.zero)
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return df_results
    end
    @error("The powerflow solver returned convergence = $(res.f_converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return res.f_converged
end

function _solve_powerflow(data::ACPowerFlowData; kwargs...)
    nlsolve_kwargs = (k for k in kwargs if first(k) ∉ AC_PF_KW)
    pf = PolarPowerFlow(data)
    J = PowerFlows.PolarPowerFlowJacobian(data, pf.x0)

    df = NLsolve.OnceDifferentiable(pf, J, pf.x0, pf.residual, J.Jv)
    res = NLsolve.nlsolve(df, pf.x0; nlsolve_kwargs...)
    return res
end
