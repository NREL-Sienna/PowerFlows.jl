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
        dc_power_flow_start!(x0, data, time_step, residual)
    else
        @debug "skipping running DC power_flow fallback"
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
function dc_power_flow_start!(x0::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    residual::ACPowerFlowResidual,
)
    _dc_power_flow_fallback!(data, time_step)
    newx0 = calculate_x0(data, time_step)
    _pick_better_x0(x0, newx0, time_step, residual, "DC power_flow fallback")
    return
end

"""Calculate x0 from data."""
function calculate_x0(data::ACPowerFlowData,
    time_step::Int64)
    n_buses = length(data.bus_type[:, 1])
    n_lcc = size(data.lcc.p_set, 1)
    x0 = Vector{Float64}(undef, 2 * n_buses + 4 * n_lcc)
    update_state!(x0, data, time_step)
    return x0
end

function _enhanced_flat_start(
    x0::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
)
    newx0 = copy(x0)
    bus_lookup = get_bus_lookup(data)
    for subnetwork_bus_axes in values(data.power_network_matrix.subnetwork_axes)
        subnetwork_indices = [bus_lookup[ix] for ix in subnetwork_bus_axes[1]]
        ref_bus = subnetwork_indices[data.bus_type[:, time_step] .== (PSY.ACBusTypes.REF,)]
        pv = subnetwork_indices[data.bus_type[:, time_step] .== (PSY.ACBusTypes.PV,)]
        pq = subnetwork_indices[data.bus_type[:, time_step] .== (PSY.ACBusTypes.PQ,)]
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
flow as a fallback. This runs a DC power_flow on `data::ACPowerFlowData` for the given
`time_step`, and writes the solution to `data.bus_angles`."""
function _dc_power_flow_fallback!(data::ACPowerFlowData, time_step::Int)
    # dev note: for DC, we can efficiently solve for all time_steps at once, and we want branch
    # flows. For AC fallback, we're only interested in the current time_step, and no branch flows
    solver_cache = get_aux_network_matrix(data).K
    # factored in constructor; no need to factor again (as long as network is same)
    valid_ix = get_valid_ix(data)
    p_inj =
        data.bus_active_power_injections[valid_ix, time_step] -
        data.bus_active_power_withdrawals[valid_ix, time_step] +
        data.bus_hvdc_net_power[valid_ix, time_step]
    # assumption: the linear algebra backend we're using implements and exports ldiv!
    ldiv!(solver_cache, p_inj)
    data.bus_angles[valid_ix, time_step] .= p_inj
end

function initialize_power_flow_variables(pf::ACPowerFlow{T},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...,
) where {T <: ACPowerFlowSolverType}
    residual = ACPowerFlowResidual(data, time_step)
    x0 = improve_x0(pf, data, residual, time_step)
    if OVERRIDE_x0 && :x0 in keys(kwargs)
        print_signorms(residual.Rv; intro = "corrected ", ps = [1, 2, Inf])
        x0 .= get(kwargs, :x0, x0)
        @warn "Overriding initial guess x0."
        residual(x0, time_step)  # re-calculate residual for new x0: might have changed.
        print_signorms(residual.Rv; ps = [1, 2, Inf])
    end
    @info "Initial residual size: " *
          "$(norm(residual.Rv, 2)) L2, " *
          "$(norm(residual.Rv, Inf)) L∞"

    J = ACPowerFlowJacobian(data, time_step)
    J(time_step)

    validate_vms::Bool = get(
        kwargs,
        :validate_voltages,
        DEFAULT_VALIDATE_VOLTAGES,
    )
    validation_range::MinMax = get(
        kwargs,
        :vm_validation_range,
        DEFAULT_VALIDATION_RANGE,
    )
    bus_types = @view get_bus_type(J.data)[:, time_step]
    validate_vms && validate_voltages(x0, bus_types, validation_range, 0)
    return residual, J, x0
end
