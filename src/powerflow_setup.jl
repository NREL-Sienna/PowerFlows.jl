function improve_x0(pf::ACPowerFlow,
    data::ACPowerFlowData,
    residual::ACPowerFlowResidual,
    time_step::Int64,
)
    x0 = calculate_x0(data, time_step)
    residual(x0, time_step)
    if norm(residual.Rv, 1) > LARGE_RESIDUAL * length(residual.Rv) &&
       get_enhanced_flat_start(pf)
        newx0 = _enhanced_flat_start(x0, data, time_step)
        _pick_better_x0(x0, newx0, time_step, residual, "enhanced flat start")
    else
        @debug "skipping enhanced flat start"
    end
    if norm(residual.Rv, 1) > LARGE_RESIDUAL * length(residual.Rv) &&
       get_robust_power_flow(pf)
        dc_powerflow_start!(x0, data, time_step, residual)
    else
        @debug "skipping running DC powerflow fallback"
    end
    residual(x0, time_step)  # re-calculate residual for new x0: might have changed.

    if sum(abs, residual.Rv) > LARGE_RESIDUAL * length(residual.Rv)
        lg_res, ix = findmax(residual.Rv)
        lg_res_rounded = round(lg_res; sigdigits = 3)
        pow_type = ix % 2 == 1 ? "active" : "reactive"
        bus_ix = div(ix + 1, 2)
        bus_no = axes(data.power_network_matrix, 1)[bus_ix]
        @warn "Initial guess provided results in a large initial residual of $lg_res_rounded. " *
              "Largest residual at bus $bus_no ($bus_ix by matrix indexing; $pow_type power)"
    end

    return x0
end

function _smaller_residual(x0::Vector{Float64},
    newx0::Vector{Float64},
    time_step::Int64,
    residual::ACPowerFlowResidual,
)
    residual(x0, time_step)
    residualSize = norm(residual.Rv, 1)
    residual(newx0, time_step)
    newResidualSize = norm(residual.Rv, 1)
    return newResidualSize < residualSize
end

function _pick_better_x0(x0::Vector{Float64},
    newx0::Vector{Float64},
    time_step::Int64,
    residual::ACPowerFlowResidual,
    improvement_method::String,
)
    if _smaller_residual(x0, newx0, time_step, residual)
        @info "success: $improvement_method yields smaller residual"
        copyto!(x0, newx0)
        residual(x0, time_step) # re-calculate for new x0.
    else
        @debug "no improvement from $improvement_method"
    end
    return nothing
end

"""If initial residual is large, run a DC power flow and see if that gives
a better starting point for angles. If so, then overwrite `x0` with the result of the DC
power flow. If not, keep the original `x0`."""
function dc_powerflow_start!(x0::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    residual::ACPowerFlowResidual,
)
    _dc_powerflow_fallback!(data, time_step)
    newx0 = calculate_x0(data, time_step)
    _pick_better_x0(x0, newx0, time_step, residual, "DC powerflow fallback")
    return
end

"""Calculate x0 from data."""
function calculate_x0(data::ACPowerFlowData,
    time_step::Int64)
    n_buses = length(data.bus_type[:, 1])
    x0 = Vector{Float64}(undef, 2 * n_buses)
    update_state!(x0, data, time_step)
    return x0
end

function _enhanced_flat_start(
    x0::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
)
    newx0 = copy(x0)
    for subnetwork_bus_axes in values(data.power_network_matrix.subnetwork_axes)
        subnetwork_indices = [data.bus_lookup[ix] for ix in subnetwork_bus_axes[1]]
        ref_bus = [
            i for
            i in subnetwork_indices if data.bus_type[i, time_step] == PSY.ACBusTypes.REF
        ]
        pv = [
            i for
            i in subnetwork_indices if data.bus_type[i, time_step] == PSY.ACBusTypes.PV
        ]
        pq = [
            i for
            i in subnetwork_indices if data.bus_type[i, time_step] == PSY.ACBusTypes.PQ
        ]
        ref_bus_angle = sum(data.bus_angles[ref_bus, time_step]) / length(ref_bus)
        if ref_bus_angle != 0.0
            newx0[2 .* vcat(pv, pq)] .= ref_bus_angle
        end
        length(pv) == 0 && length(pq) == 0 && continue
        newx0[2 .* pq .- 1] .= sum(data.bus_magnitude[pv, time_step]) / length(pv)
    end
    return newx0
end

"""When solving AC power flows, if the initial guess has large residual, we run a DC power 
flow as a fallback. This runs a DC powerflow on `data::ACPowerFlowData` for the given
`time_step`, and writes the solution to `data.bus_angles`."""
function _dc_powerflow_fallback!(data::ACPowerFlowData, time_step::Int)
    # dev note: for DC, we can efficiently solve for all timesteps at once, and we want branch
    # flows. For AC fallback, we're only interested in the current timestep, and no branch flows
    # PERF: if multi-period and multiple time steps have bad initial guesses,
    #       we're re-creating this factorization for each time step. Store it inside
    #       data.aux_network_matrix instead.
    ABA_matrix = data.aux_network_matrix.data
    solver_cache = KLULinSolveCache(ABA_matrix)
    full_factor!(solver_cache, ABA_matrix)
    p_inj =
        data.bus_activepower_injection[data.valid_ix, time_step] -
        data.bus_activepower_withdrawals[data.valid_ix, time_step]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, time_step] .= p_inj
end
