"""
    solve_powerflow!(pf::ACPowerFlow{<:ACPowerFlowSolverType}, system::PSY.System; kwargs...)

Solves the power flow in the system and writes the solution into the relevant structs.
Updates active and reactive power setpoints for generators and active and reactive
power flows for branches (calculated in the From - To direction and in the To - From direction).

Supports passing kwargs to the PF solver.

The bus types can be changed from PV to PQ if the reactive power limits are violated.

# Arguments
- `pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The power flow solver instance, can be `NewtonRaphsonACPowerFlow` or `LUACPowerFlow` (to be used for testing only).
- `system::PSY.System`: The power system model.
- `kwargs...`: Additional keyword arguments.

## Keyword Arguments
- 'check_reactive_power_limits': if `true`, the reactive power limits are enforced by changing the respective bus types from PV to PQ. Default is `false`.
- `tol`: Infinite norm of residuals under which convergence is declared. Default is `1e-9`.
- `maxIterations`: Maximum number of Newton-Raphson iterations. Default is `30`.

# Returns
- `converged::Bool`: Indicates whether the power flow solution converged.
- The power flow results are written into the system struct.

# Examples

```julia
solve_powerflow!(pf, sys)

# Passing kwargs
solve_powerflow!(pf, sys; correct_bustypes = true)

# Passing keyword arguments
solve_powerflow!(pf, sys; maxIterations=100)
```
"""
function solve_powerflow!(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    # converged must be defined in the outer scope to be visible for return
    converged = false
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            pf,
            system;
            correct_bustypes = get(kwargs, :correct_bustypes, false),
            network_reductions = get(kwargs, :network_reductions, PNM.NetworkReduction[]),
        )

        converged = solve_powerflow!(data; pf = pf, kwargs...)

        if converged
            write_powerflow_solution!(
                system,
                pf,
                data,
                get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER),
            )
            @info("PowerFlow solve converged, the results have been stored in the system")
        else
            @error("The powerflow solver returned convergence = $converged")
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
```
"""
function solve_powerflow(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    # df_results must be defined in the outer scope first to be visible for return
    df_results = Dict{String, DataFrames.DataFrame}()
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            pf,
            system;
            correct_bustypes = get(kwargs, :correct_bustypes, false),
            network_reductions = get(kwargs, :network_reductions, PNM.NetworkReduction[]),
        )

        converged = solve_powerflow!(data; pf = pf, kwargs...)

        if converged
            @info("PowerFlow solve converged, the results are exported in DataFrames")
            df_results = write_results(pf, system, data, time_step)
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
- `data::ACPowerFlowData`: The power flow data containing the grid information and initial conditions.
- `pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The power flow solver type. Defaults to `NewtonRaphsonACPowerFlow`.
- `kwargs...`: Additional keyword arguments.

# Keyword Arguments
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

    # TODO If anything in the grid topology changes, 
    #  e.g. tap positions of transformers or in service 
    #  status of branches, Yft and Ytf must be updated!
    Yft = data.power_network_matrix.branch_admittance_from_to
    Ytf = data.power_network_matrix.branch_admittance_to_from
    @assert PNM.get_bus_lookup(Yft) == get_bus_lookup(data)
    arcs = PNM.get_arc_axis(Yft)
    @assert arcs == PNM.get_arc_axis(Ytf)
    @assert length(PNM.get_bus_axis(Yft)) == length(data.bus_angles[:, 1])
    bus_lookup = get_bus_lookup(data)
    fb_ix = [bus_lookup[bus_no] for bus_no in first.(arcs)]  # from bus indices
    tb_ix = [bus_lookup[bus_no] for bus_no in last.(arcs)]   # to bus indices
    @assert length(fb_ix) == length(arcs)

    for time_step in sorted_time_steps
        converged = _ac_powerflow(data, pf, time_step; kwargs...)
        ts_converged[time_step] = converged

        if OVERWRITE_NON_CONVERGED && !converged
            # set values to NaN for not converged time steps
            data.bus_activepower_injection[:, time_step] .= NaN
            data.bus_activepower_withdrawals[:, time_step] .= NaN
            data.bus_activepower_constant_current_withdrawals[:, time_step] .= NaN
            data.bus_activepower_constant_impedance_withdrawals[:, time_step] .= NaN
            data.bus_reactivepower_injection[:, time_step] .= NaN
            data.bus_reactivepower_withdrawals[:, time_step] .= NaN
            data.bus_reactivepower_constant_current_withdrawals[:, time_step] .= NaN
            data.bus_reactivepower_constant_impedance_withdrawals[:, time_step] .= NaN
            data.bus_magnitude[:, time_step] .= NaN
            data.bus_angles[:, time_step] .= NaN
        end
    end

    # write branch flows
    # TODO if Yft, Ytf change between time steps, this must be moved inside the loop!
    # NOTE PNM's structs use ComplexF32, while the system objects store Float64's.
    #      so if you set the system bus angles/voltages to match these fields, then repeat 
    #      this math using the system voltages, you'll see differences in the flows, ~1e-4.
    ts_V =
        data.bus_magnitude[:, sorted_time_steps] .*
        exp.(1im .* data.bus_angles[:, sorted_time_steps])

    Sft = ts_V[fb_ix, :] .* conj.(Yft.data * ts_V)
    Stf = ts_V[tb_ix, :] .* conj.(Ytf.data * ts_V)
    data.arc_activepower_flow_from_to .= real.(Sft)
    data.arc_reactivepower_flow_from_to .= imag.(Sft)
    data.arc_activepower_flow_to_from .= real.(Stf)
    data.arc_reactivepower_flow_to_from .= imag.(Stf)

    data.converged .= ts_converged

    return all(data.converged)
end

function _ac_powerflow(
    data::ACPowerFlowData,
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    time_step::Int64;
    kwargs...,
)
    check_reactive_power_limits = get(kwargs, :check_reactive_power_limits, false)

    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged = _newton_powerflow(pf, data, time_step; kwargs...)
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
        Q_max = data.bus_reactivepower_bounds[ix, time_step][2]
        Q_min = data.bus_reactivepower_bounds[ix, time_step][1]

        if isnan(Q_min) && isnan(Q_max)
            @warn "Reactive power limits are uninitialized for bus $(bus_names[ix])" maxlog =
                PF_MAX_LOG
            continue
        end

        if Q_gen <= Q_min
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix, time_step] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix, time_step] =
                data.bus_reactivepower_bounds[ix, time_step][1]
        elseif Q_gen >= Q_max
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
        findall(==(bus_type), data.bus_type[:, time_step]) for bus_type in bus_types
    ]
end
