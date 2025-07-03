"""
    solve_power_flow!(system::PSY.System, ac_pf::ACPowerFlow{<:ACPowerFlowSolverType}; kwargs...)

Solves the power flow in the system and writes the solution into the relevant structs.
Updates active and reactive power setpoints for generators and active and reactive
power flows for branches (calculated in the From - To direction and in the To - From direction).

The bus types can be changed from PV to PQ if the reactive power limits are violated.

# Arguments
- `system::PSY.System`: The power system model.
- `ac_pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The power flow solver instance, can be `NewtonRaphsonACPowerFlow` or `LUACPowerFlow` (to be used for testing only).
- `kwargs...`: Additional keyword arguments.

## Keyword Arguments
- `check_connectivity::Bool`: Checks if the grid is connected. Default is `true`.
- `correct_bustypes::Bool`: If `true`, the bus types are corrected based on the generator and load data. Default is `false`.

# Returns
- `converged::Bool`: Indicates whether the power flow solution converged.
- The power flow results are written into the system struct.

# Examples

```julia
solve_power_flow!(sys, ac_pf)

# Passing kwargs
solve_power_flow!(sys, ac_pf; check_connectivity=false)
```
"""
function solve_power_flow!(
    system::PSY.System,
    ac_pf::ACPowerFlow{<:ACPowerFlowSolverType};
    kwargs...,
)
    # converged must be defined in the outer scope to be visible for return
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            ac_pf,
            system;
            check_connectivity = get(kwargs, :check_connectivity, true),
            correct_bustypes = get(kwargs, :correct_bustypes, false),
        )

        converged = _ac_powerflow(data, ac_pf, time_step)

        if converged
            write_powerflow_solution!(
                system,
                data,
                DEFAULT_NR_MAX_ITER, # could add "number of retries for reactive power limits"
                # as another kwarg.
            )
            @info("PowerFlow solve converged, the results have been stored in the system")
        else
            @error("The powerflow solver returned convergence = $(converged)")
        end
    end

    return converged
end

"""
Similar to solve_power_flow!(sys, ac_pf) but does not update the system struct with results.
Returns the results in a dictionary of dataframes.

## Examples

```julia
res = solve_power_flow(sys, ac_pf)
```
"""
function solve_power_flow(
    system::PSY.System,
    ac_pf::ACPowerFlow{<:ACPowerFlowSolverType};
    kwargs...,
)
    # df_results must be defined in the outer scope first to be visible for return
    df_results = Dict{String, DataFrames.DataFrame}()
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            ac_pf,
            system;
            check_connectivity = get(kwargs, :check_connectivity, true),
            correct_bustypes = get(kwargs, :correct_bustypes, false),
        )

        converged = _ac_powerflow(data, ac_pf, time_step)

        if converged
            @info("PowerFlow solve converged, the results are exported in DataFrames")
            df_results = write_results(ac_pf, system, data, time_step)
        else
            df_results = missing
            @error("The powerflow solver returned convergence = $(converged)")
        end
    end

    return df_results
end

"""
    solve_power_flow_data!(data::ACPowerFlowData, ac_pf::ACPowerFlow{<:ACPowerFlowSolverType})

Solve the multiperiod AC power flow problem for the given power flow data.

The bus types can be changed from PV to PQ if the reactive power limits are violated.

# Arguments
- `data::ACPowerFlowData`: The power flow data containing the grid information and initial conditions.
- `ac_pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The `ACPowerFlow` instance: this provides the solver type (Newton-Raphson, Trust Region, etc.) and its specific search parameters (tolerance, max iterations, etc.).

# Keyword Arguments
- `check_connectivity::Bool`: Checks if the grid is connected. Default is `true`.
- 'check_reactive_power_limits': if `true`, the reactive power limits are enforced by changing the respective bus types from PV to PQ. Default is `false`.
- `time_steps`: Specifies the time steps to solve. Defaults to sorting and collecting the keys of `data.timestep_map`.

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
solve_power_flow_data!(data, ACPowerFlow())
```
"""
function solve_power_flow_data!(
    data::ACPowerFlowData,
    ac_pf::ACPowerFlow{<:ACPowerFlowSolverType};
    kwargs...,
)
    # why do we need to specify time steps in 2 places like this?
    # sorted_time_steps = get(kwargs, :time_steps, sort(collect(keys(data.timestep_map))))
    sorted_time_steps = sort(collect(keys(data.timestep_map)))

    # preallocate results
    ts_converged = fill(false, length(sorted_time_steps))

    # TODO If anything in the grid topology changes, 
    #  e.g. tap positions of transformers or in service 
    #  status of branches, Yft and Ytf must be updated!
    Yft = data.power_network_matrix.yft
    Ytf = data.power_network_matrix.ytf
    fb = data.power_network_matrix.fb
    tb = data.power_network_matrix.tb

    for time_step in sorted_time_steps
        converged = _ac_powerflow(data, ac_pf, time_step)
        ts_converged[time_step] = converged

        if !converged  # set values to NaN for not converged time steps
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
    ts_V =
        data.bus_magnitude[:, sorted_time_steps] .*
        exp.(1im .* data.bus_angles[:, sorted_time_steps])
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
    ac_pf::ACPowerFlow{<:ACPowerFlowSolverType},
    time_step::Int64,
)
    check_reactive_power_limits = ac_pf.check_reactive_power_limits

    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged = _newton_powerflow(data, get_solver(ac_pf), time_step)
        if !converged || !check_reactive_power_limits ||
           _check_q_limit_bounds!(data, time_step)
            return converged
        end
    end

    @error("could not enforce reactive power limits after $MAX_REACTIVE_POWER_ITERATIONS")
    return converged
end

function _check_q_limit_bounds!(
    data::ACPowerFlowData,
    time_step::Int64,
)
    bus_names = data.power_network_matrix.axes[1]
    within_limits = true
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
        if bt == PSY.ACBusTypes.PV
            Q_gen = data.bus_reactivepower_injection[ix, time_step]
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
