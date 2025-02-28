"""
    solve_powerflow!(pf::ACPowerFlow{<:ACPowerFlowSolverType}, system::PSY.System; kwargs...)

Solves the power flow in the system and writes the solution into the relevant structs.
Updates active and reactive power setpoints for generators and active and reactive
power flows for branches (calculated in the From - To direction and in the To - From direction).

Supports passing kwargs to the PF solver.

The bus types can be changed from PV to PQ if the reactive power limits are violated.

# Arguments
- `pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The power flow solver instance, can be `NewtonRaphsonACPowerFlow` or `PowerFlows.LUACPowerFlow` (to be used for testing only).
- `system::PSY.System`: The power system model.
- `kwargs...`: Additional keyword arguments.

## Keyword Arguments
- `check_connectivity::Bool`: Checks if the grid is connected. Default is `true`.
- 'check_reactive_power_limits': if `true`, the reactive power limits are enforced by changing the respective bus types from PV to PQ. Default is `false`.
- `tol`: Infinite norm of residuals under which convergence is declared. Default is `1e-9`.
- `maxIterations`: Maximum number of Newton-Raphson iterations. Default is `30`.

# Returns
- `converged::Bool`: Indicates whether the power flow solution converged.
- The power flow results are written into the system struct.

# Examples

```julia
solve_ac_powerflow!(pf, sys)

# Passing kwargs
solve_ac_powerflow!(pf, sys; check_connectivity=false)

# Passing keyword arguments
solve_ac_powerflow!(pf, sys; maxIterations=100)
```
"""
function solve_powerflow!(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    # converged must be defined in the outer scope to be visible for return
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            pf,
            system;
            check_connectivity = get(kwargs, :check_connectivity, true),
        )

        converged, V, Sbus_result =
            _ac_powerflow(data, pf, time_step; kwargs...)
        x = _calc_x(data, V, Sbus_result, time_step)

        if converged
            write_powerflow_solution!(
                system,
                x,
                data,
                get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER),
            )
            @info("PowerFlow solve converged, the results have been stored in the system")
        else
            @error("The powerflow solver returned convergence = $(converged)")
        end
    end

    return converged
end

"""
Similar to solve_powerflow!(pf, sys) but does not update the system struct with results.
Returns the results in a dictionary of dataframes.

## Examples

```julia
res = solve_powerflow(pf, sys)

# Passing NLsolve arguments
res = solve_powerflow(pf, sys; method=:newton)
```
"""
function solve_powerflow(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    # df_results must be defined in the oueter scope first to be visible for return
    df_results = Dict{String, DataFrames.DataFrame}()
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            pf,
            system;
            check_connectivity = get(kwargs, :check_connectivity, true),
        )

        converged, V, Sbus_result = _ac_powerflow(data, pf, time_step; kwargs...)
        x = _calc_x(data, V, Sbus_result, time_step)

        if converged
            @info("PowerFlow solve converged, the results are exported in DataFrames")
            df_results = write_results(pf, system, data, x, time_step)
        else
            df_results = missing
            @error("The powerflow solver returned convergence = $(converged)")
        end
    end

    return df_results
end

"""
    solve_powerflow!(data::ACPowerFlowData; pf::ACPowerFlow{<:ACPowerFlowSolverType} = ACPowerFlow(), kwargs...)

Solve the multiperiod AC power flow problem for the given power flow data.

The bus types can be changed from PV to PQ if the reactive power limits are violated.

# Arguments
- `data::ACPowerFlowData`: The power flow data containing netwthe grid information and initial conditions.
- `pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The power flow solver type. Defaults to `NewtonRaphsonACPowerFlow`.
- `kwargs...`: Additional keyword arguments.

# Keyword Arguments
- `check_connectivity::Bool`: Checks if the grid is connected. Default is `true`.
- 'check_reactive_power_limits': if `true`, the reactive power limits are enforced by changing the respective bus types from PV to PQ. Default is `false`.
- `time_steps`: Specifies the time steps to solve. Defaults to sorting and collecting the keys of `data.timestep_map`.

# Returns
- `Nothing`. The results are written directly to the `data` object.

# Description
This function solves the AC power flow problem for each time step specified in `data`. 
It preallocates memory for the results and iterates over the sorted time steps. 
    For each time step, it calls the `_ac_powerflow` function to solve the power flow equations and updates the `data` object with the results. 
    If the power flow converges, it updates the active and reactive power injections, as well as the voltage magnitudes and angles for different bus types (REF, PV, PQ). 
    If the power flow does not converge, it sets the corresponding entries in `data` to `NaN`. 
    Finally, it calculates the branch power flows and updates the `data` object.

# Notes
- If the grid topology changes (e.g., tap positions of transformers or in-service status of branches), the admittance matrices `Yft` and `Ytf` must be updated.
- If `Yft` and `Ytf` change between time steps, the branch flow calculations must be moved inside the loop.

# Examples
```julia
solve_powerflow!(data)
```
"""
function solve_powerflow!(
    data::ACPowerFlowData;
    pf::ACPowerFlow{<:ACPowerFlowSolverType} = ACPowerFlow(),
    kwargs...,
)
    sorted_time_steps = get(kwargs, :time_steps, sort(collect(keys(data.timestep_map))))
    # preallocate results
    ts_converged = fill(false, length(sorted_time_steps))
    ts_V = zeros(Complex{Float64}, first(size(data.bus_type)), length(sorted_time_steps))
    ts_S = zeros(Complex{Float64}, first(size(data.bus_type)), length(sorted_time_steps))

    # TODO If anything in the grid topology changes, 
    #  e.g. tap positions of transformers or in service 
    #  status of branches, Yft and Ytf must be updated!
    Yft = data.power_network_matrix.yft
    Ytf = data.power_network_matrix.ytf
    fb = data.power_network_matrix.fb
    tb = data.power_network_matrix.tb

    for time_step in sorted_time_steps
        converged, V, Sbus_result =
            _ac_powerflow(data, pf, time_step; kwargs...)
        ts_converged[time_step] = converged
        ts_V[:, time_step] .= V
        ts_S[:, time_step] .= Sbus_result

        if converged
            ref, pv, pq = bus_type_idx(data, time_step)

            # write results to data:
            # write results for REF
            data.bus_activepower_injection[ref, time_step] .=
                real.(Sbus_result[ref]) .+ data.bus_activepower_withdrawals[ref, time_step]
            data.bus_reactivepower_injection[ref, time_step] .=
                imag.(Sbus_result[ref]) .+
                data.bus_reactivepower_withdrawals[ref, time_step]
            # write Q results for PV
            data.bus_reactivepower_injection[pv, time_step] .=
                imag.(Sbus_result[pv]) .+ data.bus_reactivepower_withdrawals[pv, time_step]
            # results for PQ buses do not need to be updated -> already consistent with inputs

            # write voltage results
            data.bus_magnitude[pq, time_step] .= abs.(V[pq])
            data.bus_angles[pq, time_step] .= angle.(V[pq])
            data.bus_angles[pv, time_step] .= angle.(V[pv])
        else
            data.bus_activepower_injection[:, time_step] .= NaN
            data.bus_activepower_withdrawals[:, time_step] .= NaN
            data.bus_reactivepower_injection[:, time_step] .= NaN
            data.bus_reactivepower_withdrawals[:, time_step] .= NaN
            data.bus_magnitude[:, time_step] .= NaN
            data.bus_angles[:, time_step] .= NaN
        end
    end

    # write branch flows
    # TODO if Yft, Ytf change between time steps, this must be moved inside the loop!
    Sft = ts_V[fb, :] .* conj.(Yft * ts_V)
    Stf = ts_V[tb, :] .* conj.(Ytf * ts_V)

    data.branch_activepower_flow_from_to .= real.(Sft)
    data.branch_reactivepower_flow_from_to .= imag.(Sft)
    data.branch_activepower_flow_to_from .= real.(Stf)
    data.branch_reactivepower_flow_to_from .= imag.(Stf)

    data.converged .= ts_converged

    return
end

function _ac_powerflow(
    data::ACPowerFlowData,
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    time_step::Int64;
    kwargs...,
)
    check_reactive_power_limits = get(kwargs, :check_reactive_power_limits, false)

    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged, V, Sbus_result =
            _newton_powerflow(pf, data, time_step; kwargs...)
        if !converged || !check_reactive_power_limits ||
           _check_q_limit_bounds!(data, Sbus_result, time_step)
            return (converged, V, Sbus_result)
        end
    end

    @error("could not enforce reactive power limits after $MAX_REACTIVE_POWER_ITERATIONS")
    return (converged, V, Sbus_result)
end

function _check_q_limit_bounds!(
    data::ACPowerFlowData,
    Sbus_result::Vector{Complex{Float64}},
    time_step::Int64,
)
    bus_names = data.power_network_matrix.axes[1]
    within_limits = true
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
        if bt == PSY.ACBusTypes.PV
            Q_gen = imag(Sbus_result[ix])
        else
            continue
        end

        if Q_gen <= data.bus_reactivepower_bounds[ix, time_step][1]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix, time_step] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix, time_step] =
                data.bus_reactivepower_bounds[ix, time_step][1]
        elseif Q_gen >= data.bus_reactivepower_bounds[ix, time_step][2]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix, time_step] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix, time_step] =
                data.bus_reactivepower_bounds[ix, time_step][2]
        else
            @debug "Within Limits"
        end
    end
    return within_limits
end

function _solve_powerflow!(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    data::ACPowerFlowData,
    check_reactive_power_limits,
    time_step::Int64;
    kwargs...,
)
    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged, V, Sbus_result =
            _newton_powerflow(pf, data, time_step; kwargs...)
        if !converged || !check_reactive_power_limits ||
           _check_q_limit_bounds!(data, Sbus_result, time_step)
            return converged, V, Sbus_result
        end
    end
    # todo: throw error? set converged to false? -> leave as is for now
    @error("could not enforce reactive power limits after $MAX_REACTIVE_POWER_ITERATIONS")
    return converged, V, Sbus_result
end

function bus_type_idx(
    data::ACPowerFlowData,
    time_step::Int64 = 1,
    bus_types::Tuple{Vararg{PSY.ACBusTypes}} = (
        PSY.ACBusTypes.REF,
        PSY.ACBusTypes.PV,
        PSY.ACBusTypes.PQ,
    ),
)
    # Find indices for each bus type
    return [
        findall(x -> x == bus_type, data.bus_type[:, time_step]) for bus_type in bus_types
    ]
end
