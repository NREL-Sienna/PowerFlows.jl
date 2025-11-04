function _is_available_source(x, bus::PSY.ACBus)
    # temporary workaround for FACTSControlDevice
    return PSY.get_available(x) && x.bus == bus && !isa(x, PSY.ElectricLoad) &&
           !isa(x, PSY.FACTSControlDevice)
end

"""Returns a dictionary of bus index to power contribution at that bus from FixedAdmittance
components, as a tuple of (active power, reactive power)."""
function _calculate_fixed_admittance_powers(
    sys::PSY.System,
    data::PowerFlowData,
    time_step::Int,
)
    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))
    bus_lookup = get_bus_lookup(data)

    busIxToFAPower = Dict{Int64, Tuple{Float64, Float64}}()
    for l in PSY.get_available_components(PSY.FixedAdmittance, sys)
        b = PSY.get_bus(l)
        bus_ix = PNM.get_bus_index(PSY.get_number(b), bus_lookup, nrd)
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

# sometimes errors on @assert length(remaining_unit_index) == 1. See issue #231
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

# sometimes errors on @assert length(remaining_unit_index) == 1. See issue #231
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
Calculate series voltages at buses removed in degree 2 reduction.
Method: number the nodes in the series segment 0, 1, ..., n. Number the segments by
their concluding node: 1, 2, ... n. The currents in the segments are given by:
[y^i_ff y^i_ft; y^i_tf y^i_tt] * [V_{i-1}; V_i] = [I_{i-1, i}; I_{i, i-1}]
where I'm using upper indices to denote the segment number.
There are no loads or generators at the internal nodes, so I_{i, i+1} + I_{i, i-1} = 0.
Substitute the above expressions for the currents and group by V_i:
y^i_{tf} V_{i-1} + (y_{tt}^i + y_{ff}^{i+1}) V_i + y_{ft}^{i+1} V_{i+1} = 0
For i = 1 and i = n-1, move the terms involving V_0 and V_n [known] to the other side.
This gives a tridiagonal system for x = [V_1, ..., V_{n-1}]:
A * x = [-y^1_{tf} * V_0, 0, ..., 0, -y^{n}_{ft} * V_n]
where A has diagonal entries y_{tt}^i + y_{ff}^{i+1}, subdiagonal
entries y_{tf}^{i+1}, and superdiagonal entries y_{ft}^i.

In the below, I use y_11 instead of y_ff, y_12 instead of y_ft, etc.
"""
function _set_series_voltages_and_flows!(
    sys::PSY.System,
    segment_sequence::PNM.BranchesSeries,
    equivalent_arc::Tuple{Int, Int},
    V_endpoints::Tuple{ComplexF64, ComplexF64},
    temp_bus_map::Dict{Int, String},
)
    chain_len = PNM.length(segment_sequence)
    nbuses = chain_len + 1
    # we find the voltages at interior nodes by solving Av = b, where A is tri-diagonal.
    # diagonal elements of A are: y_11 of "out" branch + y_22 of "in" branch.
    # above/below diagonal elements are y_12/y_21 of the interior segments
    d = zeros(ComplexF64, nbuses - 2)
    dl, du = zeros(ComplexF64, nbuses - 3), zeros(ComplexF64, nbuses - 3)
    b = zeros(ComplexF64, nbuses - 2)
    expected_from = equivalent_arc[1]
    y21_first, y12_last = zero(ComplexF64), zero(ComplexF64)
    for (i, segment) in enumerate(segment_sequence)
        # make sure segments are all oriented in the same direction.
        (segment_from, segment_to) = PNM.get_arc_tuple(segment)
        reversed = (segment_from != expected_from)
        @assert (!reversed) || (segment_to == expected_from)
        if !reversed
            (y11, y12, y21, y22) = PNM.ybus_branch_entries(segment)
        else
            (y11, y12, y21, y22) = reverse(PNM.ybus_branch_entries(segment))
            (segment_from, segment_to) = (segment_to, segment_from)
        end
        if i != 1 && i != chain_len
            du[i - 1] += y12
            dl[i - 1] += y21
        end
        if i != 1
            d[i - 1] += y11
        else
            y21_first = y21
        end
        if i != chain_len
            d[i] += y22
        else
            y12_last = y12
        end
        expected_from = segment_to
    end
    A = LinearAlgebra.Tridiagonal(dl, d, du)
    # if only 2 segments, these two contributions hit the same entry. thus -= c, not = -c.
    b[1] -= y21_first * V_endpoints[1]
    b[end] -= y12_last * V_endpoints[2]
    x = A \ b
    prev_bus_no, current_bus_no = equivalent_arc[1], -1
    prev_V, current_V = V_endpoints[1], zero(ComplexF64)
    # set the voltages at the interior nodes.
    # number the buses in series in order: 0, 1, 2, ... nbuses-1
    # current here is i, prev is i-1.
    for (i, segment) in enumerate(segment_sequence)
        (segment_from, segment_to) = PNM.get_arc_tuple(segment)
        reversed = segment_from != prev_bus_no
        current_bus_no = reversed ? segment_from : segment_to

        current_bus = PSY.get_component(PSY.ACBus, sys, temp_bus_map[current_bus_no])
        current_V = (i == length(segment_sequence)) ? V_endpoints[2] : x[i]
        set_voltage!(current_bus, current_V) # set voltage at bus i

        (V_from, V_to) = reversed ? (current_V, prev_V) : (prev_V, current_V)
        S = get_segment_flow(segment, V_from, V_to) # set flow at segment between i-1 and i.
        set_power_flow!(segment, S)

        prev_bus_no = current_bus_no
        if i < length(segment_sequence)
            prev_V = x[i]
        end
    end
    return
end

"""Set the power flow in the arcs that remain after network reduction. Called on the 
`direct_branch_map` and `transformer3W_map` dictionaries."""
function set_branch_flows_for_dict!(
    d::Dict{Tuple{Int, Int}, PSY.ACTransmission},
    data::ACPowerFlowData,
    time_step::Int,
)
    arc_lookup = get_arc_lookup(data)
    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))
    for (arc, br) in d
        @assert PNM.get_arc_tuple(br, nrd) == arc "disagreement between keys of " *
                                                  "map and physical arcs at arc $arc"
        @assert arc in keys(arc_lookup) "disagreement between keys of " *
                                        "map and arc axis at arc $arc"
        arc_ix = arc_lookup[arc]
        p_branch = data.arc_activepower_flow_from_to[arc_ix, time_step]
        q_branch = data.arc_reactivepower_flow_from_to[arc_ix, time_step]
        # TODO: now br could be a BranchesParallel or a ThreeWindingTransformerWinding object.
        set_power_flow!(br, p_branch + im * q_branch)
    end
end

"""
Updates system voltages and powers with power flow results
"""
function write_powerflow_solution!(
    sys::PSY.System,
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    data::ACPowerFlowData,
    max_iterations::Int,
    time_step::Int = 1,
)
    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))

    # getting bus by number is slow, O(n), so use names instead.
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )

    gspf = if isnothing(data.generator_slack_participation_factors)
        nothing
    else
        data.generator_slack_participation_factors[time_step]
    end

    # once redistribution is working again, could remove skip_redistribution.
    bus_lookup = get_bus_lookup(data)
    for (bus_number, reduced_buses) in PNM.get_bus_reduction_map(nrd)
        if length(reduced_buses) == 0
            # no reduction.
            bus_name = temp_bus_map[bus_number]
            bus = PSY.get_component(PSY.ACBus, sys, bus_name)
            ix = bus_lookup[bus_number]
            bustype = data.bus_type[ix, time_step] # may not be the same as bus.bustype!
            if bustype != PSY.get_bustype(bus)
                @warn "Changing system bus type at bus $(PSY.get_name(bus)) to match " *
                      "power flow bus type." maxlog = PF_MAX_LOG
                PSY.set_bustype!(bus, bustype)
            end
            if bustype == PSY.ACBusTypes.REF && !pf.skip_redistribution
                P_gen = data.bus_activepower_injection[ix, time_step]
                Q_gen = data.bus_reactivepower_injection[ix, time_step]
                _power_redistribution_ref(sys, P_gen, Q_gen, bus, max_iterations, gspf)
            elseif bustype == PSY.ACBusTypes.PV
                Q_gen = data.bus_reactivepower_injection[ix, time_step]
                bus.angle = data.bus_angles[ix, time_step]
                # If the PV bus has a nonzero slack participation factor,
                # then not only reactive power but also active power could have been changed
                # in the power flow calculation. This requires the same
                # active and reactive power redistribution step as for the REF bus.
                if data.bus_slack_participation_factors[ix, time_step] != 0.0 &&
                   !pf.skip_redistribution
                    P_gen = data.bus_activepower_injection[ix, time_step]
                    _power_redistribution_ref(
                        sys,
                        P_gen,
                        Q_gen,
                        bus,
                        max_iterations,
                        gspf,
                    )
                elseif !pf.skip_redistribution
                    _reactive_power_redistribution_pv(sys, Q_gen, bus, max_iterations)
                end
            elseif bustype == PSY.ACBusTypes.PQ
                Vm = data.bus_magnitude[ix, time_step]
                θ = data.bus_angles[ix, time_step]
                PSY.set_magnitude!(bus, Vm)
                PSY.set_angle!(bus, θ)
            end
        else
            @warn "Buses $reduced_buses were reduced into bus $bus_number: skipping reactive" *
                  " power redistribution and leaving voltage fields unchanged for those" *
                  " buses" maxlog = PF_MAX_LOG
        end
    end

    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))
    both_branch_types = merge(
        PNM.get_direct_branch_map(nrd),
        PNM.get_transformer3W_map(nrd),
    )
    set_branch_flows_for_dict!(
        both_branch_types,
        data,
        time_step,
    )

    if get_lcc_count(data) > 0
        # TODO LCCs and network reductions.
        arc_to_lcc = Dict{Tuple{Int, Int}, PSY.TwoTerminalLCCLine}()
        for lcc in PSY.get_available_components(PSY.TwoTerminalLCCLine, sys)
            arc_to_lcc[PNM.get_arc_tuple(PSY.get_arc(lcc))] = lcc
        end

        for (i, arc) in enumerate(data.lcc.arcs)
            lcc = arc_to_lcc[arc]
            PSY.set_rectifier_tap_setting!(lcc, data.lcc.rectifier.tap[i, time_step])
            PSY.set_inverter_tap_setting!(lcc, data.lcc.inverter.tap[i, time_step])
            PSY.set_rectifier_delay_angle!(
                lcc,
                data.lcc.rectifier.thyristor_angle[i, time_step],
            )
            PSY.set_inverter_extinction_angle!(
                lcc,
                data.lcc.inverter.thyristor_angle[i, time_step],
            )
            PSY.set_active_power_flow!(
                lcc,
                data.lcc.arc_activepower_flow_from_to[i, time_step],
            )
        end
    end

    # calculate the bus voltages at buses removed in degree 2 reduction.
    bus_lookup = get_bus_lookup(data)
    for (equivalent_arc, segments) in PNM.get_series_branch_map(nrd)
        (bus_from, bus_to) = equivalent_arc
        (ix_from, ix_to) = (bus_lookup[bus_from], bus_lookup[bus_to])
        Vm_endpoints = (data.bus_magnitude[ix_from], data.bus_magnitude[ix_to])
        Va_endpoints = (data.bus_angles[ix_from], data.bus_angles[ix_to])
        V_endpoints = Vm_endpoints .* exp.(im .* Va_endpoints)
        _set_series_voltages_and_flows!(
            sys,
            segments,
            equivalent_arc,
            V_endpoints,
            temp_bus_map,
        )
    end

    # note: this assumes all bus voltages have been written to the system objects already.
    for (equiv_arc, parallel_branches) in PNM.get_parallel_branch_map(nrd)
        (bus_from_no, bus_to_no) = equiv_arc
        (bus_from, bus_to) = (PSY.get_component(PSY.ACBus, sys, temp_bus_map[bus_from_no]),
            PSY.get_component(PSY.ACBus, sys, temp_bus_map[bus_to_no]))
        (V_from, V_to) = (bus_from, bus_to) .|> get_complex_voltage
        S = get_segment_flow(parallel_branches, V_from, V_to)
        set_power_flow!(parallel_branches, S)
    end
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
    sys_basepower::Float64,
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
    arc_activepower_losses::Vector{Float64},
    arc_reactivepower_losses::Vector{Float64},
)
    bus_df = DataFrames.DataFrame(;
        bus_number = buses,
        Vm = bus_magnitude,
        θ = bus_angles,
        P_gen = sys_basepower .* P_gen_vect,
        P_load = sys_basepower .* P_load_vect,
        P_net = sys_basepower .* (P_gen_vect - P_load_vect),
        Q_gen = sys_basepower .* Q_gen_vect,
        Q_load = sys_basepower .* Q_load_vect,
        Q_net = sys_basepower .* (Q_gen_vect - Q_load_vect),
    )
    DataFrames.sort!(bus_df, :bus_number)

    branch_df = DataFrames.DataFrame(;
        line_name = branch_names,
        bus_from = from_bus,
        bus_to = to_bus,
        P_from_to = sys_basepower .* arc_activepower_flow_from_to,
        Q_from_to = sys_basepower .* arc_reactivepower_flow_from_to,
        P_to_from = sys_basepower .* arc_activepower_flow_to_from,
        Q_to_from = sys_basepower .* arc_reactivepower_flow_to_from,
        P_losses = sys_basepower .* arc_activepower_losses,
        Q_losses = sys_basepower .* arc_reactivepower_losses,
    )
    DataFrames.sort!(branch_df, [:bus_from, :bus_to])

    return Dict("bus_results" => bus_df, "flow_results" => branch_df)
end

function add_arc_name!(arc_names::Vector{String},
    arc_names_set::Set{String},
    arc_lookup::Dict{Tuple{Int, Int}, Int},
    arc::Tuple{Int, Int},
    arc_name::String,
)
    # we don't rely on the names being unique, but it could be confusing if they aren't.
    if FORCE_UNIQUE_NAMES
        @assert !(arc_name in arc_names_set) "Arc name collision detected: $arc_name"
        push!(arc_names_set, arc_name)
    end
    arc_names[arc_lookup[arc]] = arc_name
end

"""Return the names of the arcs in the power flow data: those that correspond to branches in the system
will get the branch names, others will get a placeholder name of the form from-to."""
function get_arc_names(data::PowerFlowData)
    arc_lookup = get_arc_lookup(data)
    arc_names = fill("", length(arc_lookup))
    arc_names_set = Set(arc_names)
    # fill in names for those that directly correspond to branches in the system
    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))
    for (arc, branch) in PNM.get_direct_branch_map(nrd)
        arc_name = PSY.get_name(branch)
        add_arc_name!(arc_names, arc_names_set, arc_lookup, arc, arc_name)
    end

    # fill in transformer winding names.
    for (arc, trf_winding) in PNM.get_transformer3W_map(nrd)
        add_arc_name!(arc_names, arc_names_set, arc_lookup, arc, PNM.get_name(trf_winding))
    end
    # fill in missing names with placeholders
    for (arc, ix) in arc_lookup
        if arc_names[ix] == ""
            arc_name = "$(arc[1])-$(arc[2])"
            add_arc_name!(arc_names, arc_names_set, arc_lookup, arc, arc_name)
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
        @info "3-winding transformers included in the results export: bus-to-star flows " *
              "reported with names like 'TransformerName-primary', " *
              "'TransformerName-secondary', and 'TransformerName-tertiary'."
    end

    result_dict = Dict{Union{String, Char}, Dict{String, DataFrames.DataFrame}}()
    for i in 1:length(data.timestep_map)
        temp_dict = _allocate_results_data(
            arc_names,
            buses,
            PSY.get_base_power(sys),
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
            zeros(size(arc_names)),
            zeros(size(arc_names)),
        )
        result_dict[data.timestep_map[i]] = temp_dict
    end
    return result_dict
end

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
    busIxToFAPower = _calculate_fixed_admittance_powers(sys, data, time_step)
    for (bus_ix, fa_power) in busIxToFAPower
        data.bus_activepower_withdrawals[bus_ix, time_step] += fa_power[1]
        data.bus_reactivepower_withdrawals[bus_ix, time_step] += fa_power[2]
    end

    # NOTE: this may be different than get_bus_numbers(sys) if there's a network reduction!
    bus_numbers = PNM.get_bus_axis(data.power_network_matrix)

    arcs = PNM.get_arc_axis(data.power_network_matrix.arc_admittance_from_to)
    from_bus = first.(arcs)
    to_bus = last.(arcs)
    arc_names = get_arc_names(data)
    if length(PSY.get_components(PSY.Transformer3W, sys)) > 0
        @info "3-winding transformers included in the results export: bus-to-star flows " *
              "reported with names like 'TransformerName-primary', " *
              "'TransformerName-secondary', and 'TransformerName-tertiary'."
    end

    arc_activepower_losses =
        data.arc_activepower_flow_from_to[:, time_step] .+
        data.arc_activepower_flow_to_from[:, time_step]
    arc_reactivepower_losses =
        data.arc_reactivepower_flow_from_to[:, time_step] .+
        data.arc_reactivepower_flow_to_from[:, time_step]

    return _allocate_results_data(
        arc_names,
        bus_numbers,
        PSY.get_base_power(sys),
        from_bus,
        to_bus,
        data.bus_magnitude[:, time_step],
        data.bus_angles[:, time_step],
        data.bus_activepower_injection[:, time_step],
        data.bus_reactivepower_injection[:, time_step],
        data.bus_activepower_withdrawals[:, time_step],
        data.bus_reactivepower_withdrawals[:, time_step],
        data.arc_activepower_flow_from_to[:, time_step],
        data.arc_reactivepower_flow_from_to[:, time_step],
        data.arc_activepower_flow_to_from[:, time_step],
        data.arc_reactivepower_flow_to_from[:, time_step],
        arc_activepower_losses,
        arc_reactivepower_losses,
    )
end

"""
Modify the values in the given `System` to correspond to the given `PowerFlowData` such that
if a new `PowerFlowData` is constructed from the resulting system it is the same as `data`.
See also `write_powerflow_solution!`. NOTE that this assumes that `data` was initialized
from `sys` and then solved with no further modifications.
"""
function update_system!(sys::PSY.System, data::PowerFlowData; time_step = 1)
    nrd = PNM.get_network_reduction_data(get_power_network_matrix(data))
    if !isempty(PNM.get_reductions(nrd))
        error("update_system! does not support systems with network reductions.")
    end
    for bus in PSY.get_components(PSY.ACBus, sys)
        bus_index = get_bus_lookup(data)[PSY.get_number(bus)]
        bus_type = data.bus_type[bus_index, time_step]  # use this instead of bus.bustype to account for PV -> PQ
        if bus_type == PSY.ACBusTypes.REF
            # For REF bus, voltage and angle are fixed; update active and reactive
            P_gen = data.bus_activepower_injection[bus_index, time_step]
            Q_gen = data.bus_reactivepower_injection[bus_index, time_step]
            _power_redistribution_ref(sys, P_gen, Q_gen, bus,
                DEFAULT_MAX_REDISTRIBUTION_ITERATIONS)
        elseif bus_type == PSY.ACBusTypes.PV
            # For PV bus, active and voltage are fixed; update reactive and angle
            Q_gen = data.bus_reactivepower_injection[bus_index, time_step]
            _reactive_power_redistribution_pv(sys, Q_gen, bus,
                DEFAULT_MAX_REDISTRIBUTION_ITERATIONS)
            PSY.set_angle!(bus, data.bus_angles[bus_index, time_step])
        elseif bus_type == PSY.ACBusTypes.PQ
            # For PQ bus, active and reactive are fixed; update voltage and angle
            Vm = data.bus_magnitude[bus_index, time_step]
            PSY.set_magnitude!(bus, Vm)
            PSY.set_angle!(bus, data.bus_angles[bus_index, time_step])
            # if it used to be a PV bus, also set the Q value:
            if bus.bustype == PSY.ACBusTypes.PV
                Q_gen = data.bus_reactivepower_injection[bus_index, time_step]
                _reactive_power_redistribution_pv(sys, Q_gen, bus,
                    DEFAULT_MAX_REDISTRIBUTION_ITERATIONS)
                # now both the Q and the Vm, Va are correct for this kind of buses
            end
        end
    end
end
