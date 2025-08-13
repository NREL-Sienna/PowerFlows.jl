function _is_available_source(x, bus::PSY.ACBus)
    return PSY.get_available(x) && x.bus == bus && !isa(x, PSY.ElectricLoad)
end

"""
Calculates the From - To complex power flow (Flow injected at the bus) of branch of type
TapTransformer
"""
function flow_val(b::PSY.TapTransformer)
    !PSY.get_available(b) && return 0.0
    Y_t = PSY.get_series_admittance(b)
    c = 1 / PSY.get_tap(b)
    arc = PSY.get_arc(b)
    V_from = arc.from.magnitude * (cos(arc.from.angle) + sin(arc.from.angle) * 1im)
    V_to = arc.to.magnitude * (cos(arc.to.angle) + sin(arc.to.angle) * 1im)
    I = (V_from * Y_t * c^2) - (V_to * Y_t * c)
    flow = V_from * conj(I)
    return flow
end

"""
Calculates the From - To complex power flow (Flow injected at the bus) of branch of type
Line
"""
function flow_val(b::PSY.ACTransmission)
    !PSY.get_available(b) && return 0.0
    Y_t = PSY.get_series_admittance(b)
    arc = PSY.get_arc(b)
    V_from = arc.from.magnitude * (cos(arc.from.angle) + sin(arc.from.angle) * 1im)
    V_to = arc.to.magnitude * (cos(arc.to.angle) + sin(arc.to.angle) * 1im)
    I = V_from * (Y_t) - V_to * Y_t
    flow = V_from * conj(I)
    return flow
end

"""
Calculates the From - To complex power flow (Flow injected at the bus) of branch of type
Line
"""
function flow_val(b::PSY.DynamicBranch)
    return flow_val(b.branch)
end

"""
Calculates the From - To complex power flow (Flow injected at the bus) of branch of type
Transformer2W
"""
function flow_val(b::PSY.Transformer2W)
    !PSY.get_available(b) && return 0.0
    Y_t = PSY.get_series_admittance(b)
    arc = PSY.get_arc(b)
    V_from = arc.from.magnitude * (cos(arc.from.angle) + sin(arc.from.angle) * 1im)
    V_to = arc.to.magnitude * (cos(arc.to.angle) + sin(arc.to.angle) * 1im)
    I = V_from * (Y_t + (1im * PSY.get_primary_shunt(b))) - V_to * Y_t
    flow = V_from * conj(I)
    return flow
end

function flow_val(b::PSY.PhaseShiftingTransformer)
    error("Systems with PhaseShiftingTransformer not supported yet")
    return
end

"""
Calculates the From - To complex power flow using external data of voltages of branch of type
TapTransformer
"""
function flow_func(b::PSY.TapTransformer, V_from::Complex{Float64}, V_to::Complex{Float64})
    !PSY.get_available(b) && return (0.0, 0.0)
    Y_t = PSY.get_series_admittance(b)
    c = 1 / PSY.get_tap(b)
    I = (V_from * Y_t * c^2) - (V_to * Y_t * c)
    flow = V_from * conj(I)
    return real(flow), imag(flow)
end

"""
Calculates the From - To complex power flow using external data of voltages of branch of type
Line
"""
function flow_func(b::PSY.ACTransmission, V_from::Complex{Float64}, V_to::Complex{Float64})
    !PSY.get_available(b) && return (0.0, 0.0)
    Y_t = PSY.get_series_admittance(b)
    I = V_from * (Y_t + (1im * PSY.get_b(b).from)) - V_to * Y_t
    flow = V_from * conj(I)
    return real(flow), imag(flow)
end

"""
Calculates the From - To complex power flow using external data of voltages of branch of type
Transformer2W
"""
function flow_func(b::PSY.Transformer2W, V_from::Complex{Float64}, V_to::Complex{Float64})
    !PSY.get_available(b) && return (0.0, 0.0)
    Y_t = PSY.get_series_admittance(b)
    I = V_from * (Y_t + (1im * PSY.get_primary_shunt(b))) - V_to * Y_t
    flow = V_from * conj(I)
    return real(flow), imag(flow)
end

function flow_func(
    b::PSY.PhaseShiftingTransformer,
    V_from::Complex{Float64},
    V_to::Complex{Float64},
)
    error("Systems with PhaseShiftingTransformer not supported yet")
    return
end

"""
Updates the flow on the branches
"""
function _update_branch_flow!(sys::PSY.System)
    for b in PSY.get_components(PSY.ACTransmission, sys)
        S_flow = PSY.get_available(b) ? flow_val(b) : 0.0 + 0.0im
        PSY.set_active_power_flow!(b, real(S_flow))
        PSY.set_reactive_power_flow!(b, imag(S_flow))
    end
end

"""
Obtain total load on bus b
"""
function _get_load_data(sys::PSY.System, b::PSY.ACBus)
    active_power = 0.0
    reactive_power = 0.0
    for l in PSY.get_components(
        x -> get_available(x) && (get_bus(x) == b),
        _SingleComponentLoad,
        sys,
    )
        active_power += PSY.get_active_power(l)
        reactive_power += PSY.get_reactive_power(l)
    end
    for l in PSY.get_components(
        x -> get_available(x) && (get_bus(x) == b),
        PSY.StandardLoad,
        sys,
    )
        vm = PSY.get_magnitude(b)
        active_power +=
            PSY.get_constant_active_power(l) +
            PSY.get_current_active_power(l) * vm +
            PSY.get_impedance_active_power(l) * vm^2
        reactive_power +=
            PSY.get_constant_reactive_power(l) +
            PSY.get_current_reactive_power(l) * vm +
            PSY.get_impedance_reactive_power(l) * vm^2
    end
    return active_power, reactive_power
end

"""Returns a dictionary of bus index to power contribution at that bus from FixedAdmittance
components, as a tuple of (active power, reactive power)."""
function _calculate_fixed_admittance_powers(
    sys::PSY.System,
    data::PowerFlowData,
    time_step::Int,
)
    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))
    reverse_bus_search_map = PNM.get_reverse_bus_search_map(nrd)
    bus_lookup = get_bus_lookup(data)

    busIxToFAPower = Dict{Int64, Tuple{Float64, Float64}}()
    for l in PSY.get_components(PSY.FixedAdmittance, sys)
        !PSY.get_available(l) && continue
        b = PSY.get_bus(l)
        bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(b))
        Vm_squared =
            if get_bus_type(data)[bus_ix, time_step] == PSY.ACBusTypes.PQ
                get_bus_magnitude(data)[bus_ix, time_step]^2
            else # PV/REF bus, so V is known.
                PSY.get_magnitude(b)^2
            end
        sumSoFar = get(busIxToFAPower, bus_ix, (0.0, 0.0))
        y1, y2 = real(PSY.get_Y(l)), imag(PSY.get_Y(l))
        busIxToFAPower[bus_ix] =
            (sumSoFar[1] + y1 * Vm_squared, sumSoFar[2] - y2 * Vm_squared)
    end
    return busIxToFAPower
end

function _power_redistribution_ref(
    sys::PSY.System,
    P_gen::Float64,
    Q_gen::Float64,
    bus::PSY.ACBus,
    max_iterations::Int,
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
    } = nothing)
    devices_ =
        PSY.get_components(x -> _is_available_source(x, bus), PSY.StaticInjection, sys)
    all_devices = devices_

    sources = filter(x -> x isa PSY.Source, collect(devices_))
    non_source_devices = filter(x -> typeof(x) !== PSY.Source, collect(devices_))
    if length(sources) > 0 && length(non_source_devices) > 0
        P_gen -= sum(PSY.get_active_power.(sources))
        devices_ = setdiff(devices_, sources)
        @warn "Found sources and non-source devices at the same bus. Active power re-distribution is not well defined for this case. Source active power will remain unchanged and remaining active power will be re-distributed among non-source devices."
    elseif length(sources) > 1 && length(non_source_devices) == 0
        Psources = sum(PSY.get_active_power.(sources))
        Qsources = sum(PSY.get_reactive_power.(sources))
        if isapprox(Psources, P_gen; atol = 0.001) &&
           isapprox(Qsources, Q_gen; atol = 0.001)
            @warn "Only sources found at reference bus --- no redistribution of active or reactive power will take place"
            return
        else
            @warn "Total source P: $(Psources), Total source Q:$(Qsources) Bus P:$(P_gen), Bus Q:$(Q_gen)"
            error("Sources do not match P and/or Q requirements for reference bus.")
        end
    end
    if length(devices_) == 1
        device = first(devices_)
        PSY.set_active_power!(device, P_gen)
        _reactive_power_redistribution_pv(sys, Q_gen, bus, max_iterations)
        return
    elseif length(devices_) > 1
        devices =
            sort(collect(devices_); by = x -> get_active_power_limits_for_power_flow(x).max)
    else
        error("No devices in bus $(PSY.get_name(bus))")
    end

    if !isnothing(generator_slack_participation_factors)
        devices_gspf = Dict()
        for ((t, n), f) in generator_slack_participation_factors
            c = PSY.get_component(t, sys, n)
            PSY.get_bus(c) == bus && c ∈ all_devices && (devices_gspf[c] = f)
        end

        if isempty(devices_gspf)
            @debug "No devices with slack factors for bus $(PSY.get_name(bus))"
        else
            to_redistribute = P_gen - sum(PSY.get_active_power.(all_devices))
            sum_bus_gspf = sum(values(devices_gspf))

            for (device, factor) in devices_gspf
                PSY.set_active_power!(
                    device,
                    PSY.get_active_power(device) + to_redistribute * factor / sum_bus_gspf,
                )
            end
            _reactive_power_redistribution_pv(sys, Q_gen, bus, max_iterations)
            return
        end
    end

    sum_basepower = sum([g.max for g in get_active_power_limits_for_power_flow.(devices)])
    p_residual = P_gen
    units_at_limit = Vector{Int}()
    for (ix, d) in enumerate(devices)
        p_limits = get_active_power_limits_for_power_flow(d)
        part_factor = p_limits.max / sum_basepower
        p_frac = P_gen * part_factor
        p_set_point = clamp(p_frac, p_limits.min, p_limits.max)
        if (p_frac >= p_limits.max - BOUNDS_TOLERANCE) ||
           (p_frac <= p_limits.min + BOUNDS_TOLERANCE)
            push!(units_at_limit, ix)
            @warn "Unit $(PSY.get_name(d)) set at the limit $(p_set_point). P_max = $(p_limits.max) P_min = $(p_limits.min)"
        end
        PSY.set_active_power!(d, p_set_point)
        p_residual -= p_set_point
    end

    if !isapprox(p_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
        @debug "Ref Bus voltage residual $p_residual"
        removed_power = sum([
            g.max for g in get_active_power_limits_for_power_flow.(devices[units_at_limit])
        ])
        reallocated_p = 0.0
        it = 0
        while !isapprox(p_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
            if length(devices) == length(units_at_limit) + 1
                @warn "all devices at the active Power Limit"
                break
            end
            for (ix, d) in enumerate(devices)
                ix ∈ units_at_limit && continue
                p_limits = get_active_power_limits_for_power_flow(d)
                part_factor = p_limits.max / (sum_basepower - removed_power)
                p_frac = p_residual * part_factor
                current_p = PSY.get_active_power(d)
                p_set_point = p_frac + current_p
                if (p_set_point >= p_limits.max - BOUNDS_TOLERANCE) ||
                   (p_set_point <= p_limits.min + BOUNDS_TOLERANCE)
                    push!(units_at_limit, ix)
                    @warn "Unit $(PSY.get_name(d)) set at the limit $(p_set_point). P_max = $(p_limits.max) P_min = $(p_limits.min)"
                end
                p_set_point = clamp(p_set_point, p_limits.min, p_limits.max)
                PSY.set_active_power!(d, p_set_point)
                reallocated_p += p_frac
            end
            p_residual -= reallocated_p
            if isapprox(p_residual, 0; atol = ISAPPROX_ZERO_TOLERANCE)
                break
            end
            it += 1
            if it > max_iterations
                break
            end
        end
        if !isapprox(p_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
            remaining_unit_index = setdiff(1:length(devices), units_at_limit)
            @assert length(remaining_unit_index) == 1 remaining_unit_index
            device = devices[remaining_unit_index[1]]
            @debug "Remaining residual $q_residual, $(PSY.get_name(bus))"
            p_set_point = PSY.get_active_power(device) + p_residual
            PSY.set_active_power!(device, p_set_point)
            p_limits = get_active_power_limits_for_power_flow(device)
            if (p_set_point >= p_limits.max - BOUNDS_TOLERANCE) ||
               (p_set_point <= p_limits.min + BOUNDS_TOLERANCE)
                @error "Unit $(PSY.get_name(device)) P=$(p_set_point) above limits. P_max = $(p_limits.max) P_min = $(p_limits.min)"
            end
        end
    end
    _reactive_power_redistribution_pv(sys, Q_gen, bus, max_iterations)
    return
end

function _reactive_power_redistribution_pv(
    sys::PSY.System,
    Q_gen::Float64,
    bus::PSY.ACBus,
    max_iterations::Int,
)
    @debug "Reactive Power Distribution $(PSY.get_name(bus))"
    devices_ =
        PSY.get_components(x -> _is_available_source(x, bus), PSY.StaticInjection, sys)
    sources = filter(x -> typeof(x) == PSY.Source, collect(devices_))
    non_source_devices = filter(x -> typeof(x) !== PSY.Source, collect(devices_))
    if length(sources) > 0 && length(non_source_devices) > 0
        Q_gen -= sum(PSY.get_reactive_power.(sources))
        devices_ = setdiff(devices_, sources)
        @warn "Found sources and non-source devices at the same bus. Reactive power re-distribution is not well defined for this case. Source reactive power will remain unchanged and remaining reactive power will be re-distributed among non-source devices."
    elseif length(sources) > 1 && length(non_source_devices) == 0
        Qsources = sum(PSY.get_reactive_power.(sources))
        if isapprox(Qsources, Q_gen; atol = 0.001)
            @warn "Only sources found at PV bus --- no redistribution of reactive power will take place"
            return
        else
            @warn "Total source Q:$(Qsources), Bus Q:$(Q_gen)"
            error("Sources do not match Q requirements for PV bus.")
        end
    end
    if length(devices_) == 1
        @debug "Only one generator in the bus"
        PSY.set_reactive_power!(first(devices_), Q_gen)
        return
    elseif length(devices_) > 1
        devices = sort(collect(devices_); by = x -> PSY.get_max_reactive_power(x))
    else
        error("No devices in bus $(PSY.get_name(bus))")
    end
    # see issue https://github.com/NREL-Sienna/PowerSystems.jl/pull/1463
    total_active_power = 0.0
    for d in devices
        if PSY.get_available(d) && !isa(d, PSY.SynchronousCondenser)
            total_active_power += PSY.get_active_power(d)
        end
    end

    if isapprox(total_active_power, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
        @debug "Total Active Power Output at the bus is $(total_active_power). Using Unit's Base Power"
        sum_basepower = sum(PSY.get_base_power.(devices))
        for d in devices
            part_factor = PSY.get_base_power(d) / sum_basepower
            PSY.set_reactive_power!(d, Q_gen * part_factor)
        end
        return
    end

    q_residual = Q_gen
    units_at_limit = Vector{Int}()

    for (ix, d) in enumerate(devices)
        q_limits = get_reactive_power_limits_for_power_flow(d)
        if isapprox(q_limits.max, 0.0; atol = BOUNDS_TOLERANCE) &&
           isapprox(q_limits.min, 0.0; atol = BOUNDS_TOLERANCE)
            push!(units_at_limit, ix)
            @info "Unit $(PSY.get_name(d)) has no Q control capability. Q_max = $(q_limits.max) Q_min = $(q_limits.min)"
            continue
        end

        fraction = PSY.get_active_power(d) / total_active_power

        if fraction == 0.0
            PSY.set_reactive_power!(d, 0.0)
            continue
        else
            @assert fraction > 0
        end

        q_frac = Q_gen * fraction
        q_set_point = clamp(q_frac, q_limits.min, q_limits.max)

        if (q_frac >= q_limits.max - BOUNDS_TOLERANCE) ||
           (q_frac <= q_limits.min + BOUNDS_TOLERANCE)
            push!(units_at_limit, ix)
            @warn "Unit $(PSY.get_name(d)) set at the limit $(q_set_point). Q_max = $(q_limits.max) Q_min = $(q_limits.min)"
        end

        PSY.set_reactive_power!(d, q_set_point)
        q_residual -= q_set_point

        if isapprox(q_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
            break
        end
    end

    if !isapprox(q_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
        it = 0
        while !isapprox(q_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
            if length(devices) == length(units_at_limit) + 1
                @debug "Only one device not at the limit in Bus"
                break
            end
            removed_power = sum(PSY.get_active_power.(devices[units_at_limit]))
            reallocated_q = 0.0
            for (ix, d) in enumerate(devices)
                ix ∈ units_at_limit && continue
                q_limits = get_reactive_power_limits_for_power_flow(d)

                if removed_power < total_active_power
                    fraction =
                        PSY.get_active_power(d) / (total_active_power - removed_power)
                elseif isapprox(removed_power, total_active_power)
                    fraction = 1
                else
                    error("Remove power can't be larger than the total active power")
                end

                if fraction == 0.0
                    continue
                else
                    PSY.InfrastructureSystems.@assert_op fraction > 0
                end
                current_q = PSY.get_reactive_power(d)
                q_frac = q_residual * fraction
                q_set_point = clamp(q_frac + current_q, q_limits.min, q_limits.max)
                # Assign new capacity based on the limits and the fraction
                reallocated_q += q_set_point - current_q
                if ((q_frac + current_q) >= q_limits.max - BOUNDS_TOLERANCE) ||
                   ((q_frac + current_q) <= q_limits.min + BOUNDS_TOLERANCE)
                    push!(units_at_limit, ix)
                    @warn "Unit $(PSY.get_name(d)) set at the limit $(q_set_point). Q_max = $(q_limits.max) Q_min = $(q_limits.min)"
                end

                PSY.set_reactive_power!(d, q_set_point)
            end
            q_residual -= reallocated_q
            if isapprox(q_residual, 0; atol = ISAPPROX_ZERO_TOLERANCE)
                break
            end
            it += 1
            if it > max_iterations
                @warn "Maximum number of iterations for Q-redistribution reached. Number of devices at Q limit are: $(length(units_at_limit)) of $(length(devices)) available devices"
                break
            end
        end
    end

    # Last attempt to allocate reactive power
    if !isapprox(q_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
        remaining_unit_index = setdiff(1:length(devices), units_at_limit)
        @assert length(remaining_unit_index) == 1 remaining_unit_index
        device = devices[remaining_unit_index[1]]
        @debug "Remaining residual $q_residual, $(PSY.get_name(bus))"
        q_set_point = PSY.get_reactive_power(device) + q_residual
        PSY.set_reactive_power!(device, q_set_point)
        q_limits = get_reactive_power_limits_for_power_flow(device)
        if (q_set_point >= q_limits.max - BOUNDS_TOLERANCE) ||
           (q_set_point <= q_limits.min + BOUNDS_TOLERANCE)
            @error "Unit $(PSY.get_name(device)) Q=$(q_set_point) above limits. Q_max = $(q_limits.max) Q_min = $(q_limits.min)"
        end
    end

    @assert isapprox(
        sum(PSY.get_reactive_power.(devices)),
        Q_gen;
        atol = ISAPPROX_ZERO_TOLERANCE,
    )

    return
end

"""
Updates system voltages and powers with power flow results
"""
function write_powerflow_solution!(
    sys::PSY.System,
    data::PowerFlowData,
    max_iterations::Int,
    time_step::Int = 1,
)
    buses = enumerate(
        sort!(collect(PSY.get_components(PSY.ACBus, sys)); by = x -> PSY.get_number(x)),
    )

    # Handle any changes made manually to the PowerFlowData, not necessarily reflected in the solver result
    # Right now the only such change we handle is the one in _check_q_limit_bounds!
    for (ix, bus) in buses
        system_bustype = PSY.get_bustype(bus)
        data_bustype = data.bus_type[ix, time_step]
        (system_bustype == data_bustype) && continue
        @assert system_bustype == PSY.ACBusTypes.PV
        @assert data_bustype == PSY.ACBusTypes.PQ
        Q_gen = data.bus_reactivepower_injection[ix, time_step]
        @debug "Updating bus $(PSY.get_name(bus)) reactive power and type to PQ due to check_reactive_power_limits: $Q_gen"
        _reactive_power_redistribution_pv(sys, Q_gen, bus, max_iterations)
        PSY.set_bustype!(bus, data_bustype)
    end

    gspf = if isnothing(data.generator_slack_participation_factors)
        nothing
    else
        data.generator_slack_participation_factors[time_step]
    end

    for (ix, bus) in buses
        if bus.bustype == PSY.ACBusTypes.REF
            P_gen = data.bus_activepower_injection[ix, time_step]
            Q_gen = data.bus_reactivepower_injection[ix, time_step]
            _power_redistribution_ref(sys, P_gen, Q_gen, bus, max_iterations, gspf)
        elseif bus.bustype == PSY.ACBusTypes.PV
            Q_gen = data.bus_reactivepower_injection[ix, time_step]
            bus.angle = data.bus_angles[ix, time_step]
            # If the PV bus has a nonzero slack participation factor,
            # then not only reactive power but also active power could have been changed
            # in the power flow calculation. This requires the same
            # active and reactive power redistribution step as for the REF bus.
            if data.bus_slack_participation_factors[ix, time_step] != 0.0
                P_gen = data.bus_activepower_injection[ix, time_step]
                _power_redistribution_ref(sys, P_gen, Q_gen, bus, max_iterations, gspf)
            else
                _reactive_power_redistribution_pv(sys, Q_gen, bus, max_iterations)
            end
        elseif bus.bustype == PSY.ACBusTypes.PQ
            Vm = data.bus_magnitude[ix, time_step]
            θ = data.bus_angles[ix, time_step]
            PSY.set_magnitude!(bus, Vm)
            PSY.set_angle!(bus, θ)
        end
    end

    _update_branch_flow!(sys)
    return
end

# returns list of arc tuples and buses numbers: ABA case
function _get_arcs_buses(data::ABAPowerFlowData)
    return axes(data.aux_network_matrix)[2], axes(data.aux_network_matrix)[1]
end

# returns list of arc tuples and buses numbers: PTDF and virtual PTDF case
function _get_arcs_buses(data::Union{PTDFPowerFlowData, vPTDFPowerFlowData})
    return PNM.get_arc_axis(data.power_network_matrix),
    PNM.get_bus_axis(data.power_network_matrix)
end

function _allocate_results_data(
    branch_names::Vector{String},
    buses::Vector{Int64},
    from_bus::Vector{Int64},
    to_bus::Vector{Int64},
    bus_magnitude::Vector{Float64},
    bus_angles::Vector{Float64},
    P_gen_vect::Vector{Float64},
    Q_gen_vect::Vector{Float64},
    P_load_vect::Vector{Float64},
    Q_load_vect::Vector{Float64},
    arc_activepower_flow_from_to::Vector{Float64},
    arc_reactivepower_flow_from_to::Vector{Float64},
    arc_activepower_flow_to_from::Vector{Float64},
    arc_reactivepower_flow_to_from::Vector{Float64},
)
    bus_df = DataFrames.DataFrame(;
        bus_number = buses,
        Vm = bus_magnitude,
        θ = bus_angles,
        P_gen = P_gen_vect,
        P_load = P_load_vect,
        P_net = P_gen_vect - P_load_vect,
        Q_gen = Q_gen_vect,
        Q_load = Q_load_vect,
        Q_net = Q_gen_vect - Q_load_vect,
    )
    DataFrames.sort!(bus_df, :bus_number)

    branch_df = DataFrames.DataFrame(;
        line_name = branch_names,
        bus_from = from_bus,
        bus_to = to_bus,
        P_from_to = arc_activepower_flow_from_to,
        Q_from_to = arc_reactivepower_flow_from_to,
        P_to_from = arc_activepower_flow_to_from,
        Q_to_from = arc_reactivepower_flow_to_from,
        P_losses = zeros(length(branch_names)),
        Q_losses = zeros(length(branch_names)),
    )
    DataFrames.sort!(branch_df, [:bus_from, :bus_to])

    return Dict("bus_results" => bus_df, "flow_results" => branch_df)
end

"""Return the names of the arcs in the power flow data: those that correspond to branches in the system
will get the branch names, others will get a placeholder name of the form from-to."""
function get_arc_names(data::PowerFlowData)
    arc_lookup = get_arc_lookup(data)
    arc_names = fill("", length(arc_lookup))
    # fill in names for those that directly correspond to branches in the system
    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))
    for (arc, branch) in PNM.get_direct_branch_map(nrd)
        arc_names[arc_lookup[arc]] = PSY.get_name(branch)
    end
    # fill in missing names with placeholders
    for (arc, ix) in arc_lookup
        if arc_names[ix] == ""
            arc_names[ix] = "$(arc[1])-$(arc[2])"
        end
    end
    return arc_names
end

"""
Returns a dictionary containing the DC power flow results. Each key corresponds
to the name of the considered time periods, storing a DataFrame with the PF
results.

# Arguments:
- `data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData}`:
        PowerFlowData structure containing power flows and bus angles.
- `sys::PSY.System`:
        container storing the system information.
"""
function write_results(
    data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData},
    sys::PSY.System,
)
    @info("Voltages are exported in pu. Powers are exported in MW/MVAr.")

    ### non time-dependent variables

    # get bus and arcs
    arcs, buses = _get_arcs_buses(data)

    from_bus = first.(arcs)
    to_bus = last.(arcs)
    arc_names = get_arc_names(data)
    if length(PSY.get_components(PSY.Transformer3W, sys)) > 0
        @warn "3-winding transformers are not supported in the results export"
    end

    result_dict = Dict{Union{String, Char}, Dict{String, DataFrames.DataFrame}}()
    for i in 1:length(data.timestep_map)
        temp_dict = _allocate_results_data(
            arc_names,
            buses,
            from_bus,
            to_bus,
            data.bus_magnitude[:, i],
            data.bus_angles[:, i],
            data.bus_activepower_injection[:, i],
            data.bus_reactivepower_injection[:, i],
            data.bus_activepower_withdrawals[:, i],
            data.bus_reactivepower_withdrawals[:, i],
            data.arc_activepower_flow_from_to[:, i],
            data.arc_reactivepower_flow_from_to[:, i],
            data.arc_activepower_flow_to_from[:, i],
            data.arc_reactivepower_flow_to_from[:, i],
        )
        result_dict[data.timestep_map[i]] = temp_dict
    end
    return result_dict
end

# TODO: multi-period still to implement
"""
Returns a dictionary containing the AC power flow results.

Only single-period evaluation is supported at the moment for AC Power flows. Resulting
dictionary will therefore feature just one key linked to one DataFrame.

# Arguments:
- `::ACPowerFlow`:
        use ACPowerFlow() storing AC power flow results.
- `sys::PSY.System`:
        container storing the system information.
- `result::Vector{Float64}`:
        vector containing the results for one single time-period.
"""
function write_results(
    ::ACPowerFlow{<:ACPowerFlowSolverType},
    sys::PSY.System,
    data::ACPowerFlowData,
    time_step::Int64,
)
    @info("Voltages are exported in pu. Powers are exported in MW/MVAr.")
    sys_basepower = PSY.get_base_power(sys)
    Vm_vect = data.bus_magnitude[:, time_step]
    θ_vect = data.bus_angles[:, time_step]
    P_gen_vect = sys_basepower .* data.bus_activepower_injection[:, time_step]
    Q_gen_vect = sys_basepower .* data.bus_reactivepower_injection[:, time_step]
    P_load_vect = sys_basepower .* data.bus_activepower_withdrawals[:, time_step]
    Q_load_vect = sys_basepower .* data.bus_reactivepower_withdrawals[:, time_step]

    busIxToFAPower = _calculate_fixed_admittance_powers(sys, data, time_step)
    for (bus_ix, fa_power) in busIxToFAPower
        P_load_vect[bus_ix] += fa_power[1] * sys_basepower
        Q_load_vect[bus_ix] += fa_power[2] * sys_basepower
    end

    # NOTE: this may be different than get_bus_numbers(sys) if there's a network reduction!
    bus_numbers = PNM.get_bus_axis(data.power_network_matrix)

    bus_df = DataFrames.DataFrame(;
        bus_number = bus_numbers,
        Vm = Vm_vect,
        θ = θ_vect,
        P_gen = P_gen_vect,
        P_load = P_load_vect,
        P_net = P_gen_vect - P_load_vect,
        Q_gen = Q_gen_vect,
        Q_load = Q_load_vect,
        Q_net = Q_gen_vect - Q_load_vect,
    )

    # TODO Ybus lacks a get_arc_axis function.
    # it really should have one: the rows of the yft/ytf arrays correspond to arcs.
    arcs = fill((0, 0), length(get_arc_lookup(data)))
    for (arc, ix) in get_arc_lookup(data)
        arcs[ix] = arc
    end

    from_bus = first.(arcs)
    to_bus = last.(arcs)
    arc_names = get_arc_names(data)
    if length(PSY.get_components(PSY.Transformer3W, sys)) > 0
        @warn "3-winding transformers are not supported in the results export"
    end

    branch_df = DataFrames.DataFrame(;
        line_name = arc_names,
        bus_from = from_bus,
        bus_to = to_bus,
        P_from_to = data.arc_activepower_flow_from_to[:, time_step],
        Q_from_to = data.arc_reactivepower_flow_from_to[:, time_step],
        P_to_from = data.arc_activepower_flow_to_from[:, time_step],
        Q_to_from = data.arc_reactivepower_flow_to_from[:, time_step],
        P_losses = data.arc_activepower_flow_from_to[:, time_step] .+
                   data.arc_activepower_flow_to_from[:, time_step],
        Q_losses = data.arc_reactivepower_flow_from_to[:, time_step] .+
                   data.arc_reactivepower_flow_to_from[:, time_step],
    )
    DataFrames.sort!(branch_df, [:bus_from, :bus_to])

    return Dict("bus_results" => bus_df, "flow_results" => branch_df)
end

"""
Modify the values in the given `System` to correspond to the given `PowerFlowData` such that
if a new `PowerFlowData` is constructed from the resulting system it is the same as `data`.
See also `write_powerflow_solution!`. NOTE that this assumes that `data` was initialized
from `sys` and then solved with no further modifications.
"""
function update_system!(sys::PSY.System, data::PowerFlowData; time_step = 1)
    # TODO: throw an error if there's a network reduction.
    for bus in PSY.get_components(PSY.ACBus, sys)
        bus_number = data.bus_lookup[PSY.get_number(bus)]
        bus_type = data.bus_type[bus_number, time_step]  # use this instead of bus.bustype to account for PV -> PQ
        if bus_type == PSY.ACBusTypes.REF
            # For REF bus, voltage and angle are fixed; update active and reactive
            P_gen = data.bus_activepower_injection[bus_number, time_step]
            Q_gen = data.bus_reactivepower_injection[bus_number, time_step]
            _power_redistribution_ref(sys, P_gen, Q_gen, bus,
                DEFAULT_MAX_REDISTRIBUTION_ITERATIONS)
        elseif bus_type == PSY.ACBusTypes.PV
            # For PV bus, active and voltage are fixed; update reactive and angle
            Q_gen = data.bus_reactivepower_injection[bus_number, time_step]
            _reactive_power_redistribution_pv(sys, Q_gen, bus,
                DEFAULT_MAX_REDISTRIBUTION_ITERATIONS)
            PSY.set_angle!(bus, data.bus_angles[bus_number, time_step])
        elseif bus_type == PSY.ACBusTypes.PQ
            # For PQ bus, active and reactive are fixed; update voltage and angle
            Vm = data.bus_magnitude[bus_number, time_step]
            PSY.set_magnitude!(bus, Vm)
            PSY.set_angle!(bus, data.bus_angles[bus_number, time_step])
            # if it used to be a PV bus, also set the Q value:
            if bus.bustype == PSY.ACBusTypes.PV
                Q_gen = data.bus_reactivepower_injection[bus_number, time_step]
                _reactive_power_redistribution_pv(sys, Q_gen, bus,
                    DEFAULT_MAX_REDISTRIBUTION_ITERATIONS)
                # now both the Q and the Vm, Va are correct for this kind of buses
            end
        end
    end
end
