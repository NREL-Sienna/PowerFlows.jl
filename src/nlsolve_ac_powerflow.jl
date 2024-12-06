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
solve_ac_powerflow!(sys)
# Passing NLsolve arguments
solve_ac_powerflow!(sys, method=:newton)
```
"""
function solve_ac_powerflow!(pf::Union{NLSolveACPowerFlow,KLUACPowerFlow}, system::PSY.System; kwargs...)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    check_reactive_power_limits = get(kwargs, :check_reactive_power_limits, false)
    data = PowerFlowData(
        NLSolveACPowerFlow(; check_reactive_power_limits=check_reactive_power_limits),
        system;
        check_connectivity=get(kwargs, :check_connectivity, true),
    )
    max_iterations = DEFAULT_MAX_REDISTRIBUTION_ITERATIONS
    converged, x = _solve_powerflow!(pf, data, check_reactive_power_limits; kwargs...)
    if converged
        write_powerflow_solution!(system, x, max_iterations)
        @info("PowerFlow solve converged, the results have been stored in the system")
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return converged
    end
    @error("The powerflow solver returned convergence = $(converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return converged
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
    pf::Union{NLSolveACPowerFlow,KLUACPowerFlow},
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
        check_connectivity=get(kwargs, :check_connectivity, true),
    )

    converged, x = _solve_powerflow!(pf, data, pf.check_reactive_power_limits; kwargs...)

    if converged
        @info("PowerFlow solve converged, the results are exported in DataFrames")
        df_results = write_results(pf, system, data, x)
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return df_results
    end
    @error("The powerflow solver returned convergence = $(converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return converged
end

function _check_q_limit_bounds!(data::ACPowerFlowData, zero::Vector{Float64})
    bus_names = data.power_network_matrix.axes[1]
    within_limits = true
    for (ix, b) in enumerate(data.bus_type)
        if b == PSY.ACBusTypes.PV
            Q_gen = zero[2*ix-1]
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
    pf::Union{NLSolveACPowerFlow,KLUACPowerFlow},
    data::ACPowerFlowData,
    check_reactive_power_limits;
    nlsolve_kwargs...,
)
    if check_reactive_power_limits
        for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
            converged, x = _nlsolve_powerflow(pf, data; nlsolve_kwargs...)
            if converged
                if _check_q_limit_bounds!(data, x)
                    return converged, x
                end
            else
                return converged, x
            end
        end
    else
        return _nlsolve_powerflow(pf, data; nlsolve_kwargs...)
    end
end

function _nlsolve_powerflow(pf::NLSolveACPowerFlow, data::ACPowerFlowData; nlsolve_kwargs...)
    pf = PolarPowerFlow(data)
    J = PowerFlows.PolarPowerFlowJacobian(data, pf.x0)

    df = NLsolve.OnceDifferentiable(pf, J, pf.x0, pf.residual, J.Jv)
    res = NLsolve.nlsolve(df, pf.x0; nlsolve_kwargs...)
    if !res.f_converged
        @error("The powerflow solver returned convergence = $(res.f_converged)")
    end
    return res.f_converged, res.zero
end

function _nlsolve_powerflow(pf::KLUACPowerFlow, data::ACPowerFlowData; nlsolve_kwargs...)
    pf = PolarPowerFlow(data)
    #J_function = PowerFlows.PolarPowerFlowJacobian(data, pf.x0)

    maxIter = 30  # TODO
    tol = 1e-6    # TODO
    i = 0

    Vm = data.bus_magnitude[:]
    Va = data.bus_angles[:]
    V = Vm .* exp.(1im * Va)

    # Find indices for each bus type
    ref = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.REF, data.bus_type)
    pv = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PV, data.bus_type)
    pq = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PQ, data.bus_type)

    Ybus = pf.data.power_network_matrix.data

    Sbus = data.bus_activepower_injection[:] - data.bus_activepower_withdrawals[:] + 1im * (data.bus_reactivepower_injection[:] - data.bus_reactivepower_withdrawals[:])

    mis = V .* conj(Ybus * V) - Sbus
    F = [real(mis[[pv; pq]]); imag(mis[pq])]

    # nref = length(ref)
    npv = length(pv)
    npq = length(pq)

    converged = (npv + npq) == 0  # if only ref buses present, we do not need to enter the loop

    while i < maxIter && !converged
        i += 1
        diagV = LinearAlgebra.Diagonal(V)
        diagIbus = LinearAlgebra.Diagonal(Ybus * V)
        diagVnorm = LinearAlgebra.Diagonal(V ./ abs.(V))
        dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
        dSbus_dVa = 1im * diagV * conj(diagIbus - Ybus * diagV)

        j11 = real(dSbus_dVa[[pv; pq], [pv; pq]])
        j12 = real(dSbus_dVm[[pv; pq], pq])
        j21 = imag(dSbus_dVa[pq, [pv; pq]])
        j22 = imag(dSbus_dVm[pq, pq])
        J = [j11 j12; j21 j22]

        factor_J = KLU.klu(J)
        dx = -(factor_J \ F)

        #J_function(J_function.Jv, x)
        #dx = - (J_function.Jv \ F)

        Va[pv] .+= dx[1:npv]
        Va[pq] .+= dx[(npv+1):(npv+npq)]
        Vm[pq] .+= dx[(npv+npq+1):(npv+2*npq)]

        V = Vm .* exp.(1im * Va)

        Vm = abs.(V)
        Va = angle.(V)

        mis = V .* conj(Ybus * V) - Sbus
        F = [real(mis[[pv; pq]]); imag(mis[pq])]
        converged = LinearAlgebra.norm(F, Inf) < tol
    end


    if !converged
        @error("The powerflow solver with KLU did not converge after $i iterations")
    else
        @debug("The powerflow solver with KLU converged after $i iterations")
    end

    # mock the expected x format, where the values depend on the type of the bus:
    n_buses = length(data.bus_type)
    x = Float64[0.0 for _ in 1:(2*n_buses)]
    Sbus_result = V .* conj(Ybus * V)
    for (ix, b) in enumerate(data.bus_type)
        if b == PSY.ACBusTypes.REF
            # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
            x[2*ix-1] = real(Sbus_result[ix]) + data.bus_activepower_withdrawals[ix]
            x[2*ix] = imag(Sbus_result[ix]) + data.bus_reactivepower_withdrawals[ix]
        elseif b == PSY.ACBusTypes.PV
            # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
            x[2*ix-1] = imag(Sbus_result[ix]) + data.bus_reactivepower_withdrawals[ix]
            x[2*ix] = Va[ix]
        elseif b == PSY.ACBusTypes.PQ
            # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
            x[2*ix-1] = Vm[ix]
            x[2*ix] = Va[ix]
        end
    end

    return converged, x
end
