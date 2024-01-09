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
    check_reactive_power_limits = get(kwargs, :check_reactive_power_limits, false)
    data = PowerFlowData(
        ACPowerFlow(; check_reactive_power_limits = check_reactive_power_limits),
        system;
        check_connectivity = get(kwargs, :check_connectivity, true),
    )
    res = _solve_powerflow!(data, check_reactive_power_limits; kwargs...)
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
    pf::ACPowerFlow,
    system::PSY.System;
    kwargs...,
)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    data = PowerFlowData(
        pf,
        system;
        check_connectivity = get(kwargs, :check_connectivity, true),
    )

    res = _solve_powerflow!(data, pf.check_reactive_power_limits; kwargs...)

    if res.f_converged
        @info("PowerFlow solve converged, the results are exported in DataFrames")
        df_results = write_results(pf, system, data, res.zero)
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return df_results
    end
    @error("The powerflow solver returned convergence = $(res.f_converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return res.f_converged
end

function _check_q_limit_bounds!(data::ACPowerFlowData, zero::Vector{Float64})
    bus_names = data.power_network_matrix.axes[1]
    within_limits = true
    for (ix, b) in enumerate(data.bus_type)
        if b == PSY.ACBusTypes.PV
            Q_gen = zero[2 * ix - 1]
        else
            continue
        end

        if Q_gen <= data.bus_reactivepower_bounds[ix][1]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix] = data.bus_reactivepower_bounds[ix][1]
        elseif Q_gen >= data.bus_reactivepower_bounds[ix][2]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix] = data.bus_reactivepower_bounds[ix][2]
        else
            @debug "Within Limits"
        end
    end
    return within_limits
end

function _solve_powerflow!(
    data::ACPowerFlowData,
    check_reactive_power_limits;
    nlsolve_kwargs...,
)
    if check_reactive_power_limits
        for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
            res = _nlsolve_powerflow(data; nlsolve_kwargs...)
            if res.f_converged
                if _check_q_limit_bounds!(data, res.zero)
                    return res
                end
            else
                return res
            end
        end
    else
        return _nlsolve_powerflow(data; nlsolve_kwargs...)
    end
end

function _nlsolve_powerflow(data::ACPowerFlowData; nlsolve_kwargs...)
    pf = PolarPowerFlow(data)
    J = PolarPowerFlowJacobian(data, pf.x0)

    df = NLsolve.OnceDifferentiable(pf, J, pf.x0, pf.residual, J.Jv)
    res = NLsolve.nlsolve(df, pf.x0; nlsolve_kwargs...)
    if !res.f_converged
        @error("The powerflow solver returned convergence = $(res.f_converged)")
    end
    return res
end
