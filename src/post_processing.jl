function _is_available_source(x, bus::PSY.Bus)
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
function flow_val(b::PSY.ACBranch)
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
function flow_func(b::PSY.ACBranch, V_from::Complex{Float64}, V_to::Complex{Float64})
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
    for b in PSY.get_components(PSY.ACBranch, sys)
        S_flow = PSY.get_available(b) ? flow_val(b) : 0.0 + 0.0im
        PSY.set_active_power_flow!(b, real(S_flow))
        PSY.set_reactive_power_flow!(b, imag(S_flow))
    end
end

"""
Obtain total load on bus b
"""
function _get_load_data(sys::PSY.System, b::PSY.Bus)
    active_power = 0.0
    reactive_power = 0.0
    for l in PSY.get_components(x -> !isa(x, PSY.FixedAdmittance), PSY.ElectricLoad, sys)
        !PSY.get_available(l) && continue
        if (l.bus == b)
            active_power += get_total_p(l)
            reactive_power += get_total_q(l)
        end
    end
    return active_power, reactive_power
end

function _get_fixed_admittance_power(
    sys::PSY.System,
    b::PSY.Bus,
    result::AbstractVector,
    ix::Int,
)
    active_power = 0.0
    reactive_power = 0.0
    for l in PSY.get_components(PSY.FixedAdmittance, sys)
        !PSY.get_available(l) && continue
        if (l.bus == b)
            Vm_squared =
                b.bustype == PSY.BusTypes.PQ ? result[2 * ix - 1]^2 : PSY.get_magnitude(b)^2
            active_power += Vm_squared * real(PSY.get_Y(l))
            reactive_power -= Vm_squared * imag(PSY.get_Y(l))
        end
    end
    return active_power, reactive_power
end

function _get_limits_for_power_distribution(gen::PSY.StaticInjection)
    return PSY.get_active_power_limits(gen)
end

function _get_limits_for_power_distribution(gen::PSY.RenewableDispatch)
    return (min = 0.0, max = PSY.get_max_active_power(gen))
end

function _get_limits_for_power_distribution(gen::PSY.Storage)
    return (min = 0.0, max = PSY.get_output_active_power_limits(gen).max)
end

function _power_redistribution_ref(
    sys::PSY.System,
    P_gen::Float64,
    Q_gen::Float64,
    bus::PSY.Bus,
)
    devices_ =
        PSY.get_components(x -> _is_available_source(x, bus), PSY.StaticInjection, sys)
    sources = filter(x -> typeof(x) == PSY.Source, collect(devices_))
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
        _reactive_power_redistribution_pv(sys, Q_gen, bus)
        return
    elseif length(devices_) > 1
        devices =
            sort(collect(devices_); by = x -> _get_limits_for_power_distribution(x).max)
    else
        error("No devices in bus $(PSY.get_name(bus))")
    end

    sum_basepower = sum([g.max for g in _get_limits_for_power_distribution.(devices)])
    p_residual = P_gen
    units_at_limit = Vector{Int}()
    for (ix, d) in enumerate(devices)
        p_limits = _get_limits_for_power_distribution(d)
        part_factor = p_limits.max / sum_basepower
        p_frac = P_gen * part_factor
        p_set_point = clamp(p_frac, p_limits.min, p_limits.max)
        if (p_frac >= p_limits.max + BOUNDS_TOLERANCE) ||
           (p_frac <= p_limits.min - BOUNDS_TOLERANCE)
            push!(units_at_limit, ix)
            @warn "Unit $(PSY.get_name(d)) set at the limit $(p_set_point). P_max = $(p_limits.max) P_min = $(p_limits.min)"
        end
        PSY.set_active_power!(d, p_set_point)
        p_residual -= p_set_point
    end

    if !isapprox(p_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
        @debug "Ref Bus voltage residual $p_residual"
        removed_power = sum([
            g.max for g in _get_limits_for_power_distribution.(devices[units_at_limit])
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
                p_limits = PSY.get_active_power_limits(d)
                part_factor = p_limits.max / (sum_basepower - removed_power)
                p_frac = p_residual * part_factor
                current_p = PSY.get_active_power(d)
                p_set_point = p_frac + current_p
                if (p_set_point >= p_limits.max + BOUNDS_TOLERANCE) ||
                   (p_set_point <= p_limits.min - BOUNDS_TOLERANCE)
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
            if it > 10
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
            p_limits = PSY.get_reactive_power_limits(device)
            if (p_set_point >= p_limits.max + BOUNDS_TOLERANCE) ||
               (p_set_point <= p_limits.min - BOUNDS_TOLERANCE)
                @error "Unit $(PSY.get_name(device)) P=$(p_set_point) above limits. P_max = $(p_limits.max) P_min = $(p_limits.min)"
            end
        end
    end
    _reactive_power_redistribution_pv(sys, Q_gen, bus)
    return
end

function _reactive_power_redistribution_pv(sys::PSY.System, Q_gen::Float64, bus::PSY.Bus)
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

    total_active_power = sum(PSY.get_active_power.(devices))

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
        q_limits = PSY.get_reactive_power_limits(d)
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

        if (q_frac >= q_limits.max + BOUNDS_TOLERANCE) ||
           (q_frac <= q_limits.min - BOUNDS_TOLERANCE)
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
                q_limits = PSY.get_reactive_power_limits(d)

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
                reallocated_q += q_frac
                if (q_frac >= q_limits.max + BOUNDS_TOLERANCE) ||
                   (q_frac <= q_limits.min - BOUNDS_TOLERANCE)
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
            if it > 1
                break
            end
        end
    end

    if !isapprox(q_residual, 0.0; atol = ISAPPROX_ZERO_TOLERANCE)
        remaining_unit_index = setdiff(1:length(devices), units_at_limit)
        @assert length(remaining_unit_index) == 1 remaining_unit_index
        device = devices[remaining_unit_index[1]]
        @debug "Remaining residual $q_residual, $(PSY.get_name(bus))"
        q_set_point = PSY.get_reactive_power(device) + q_residual
        PSY.set_reactive_power!(device, q_set_point)
        q_limits = PSY.get_reactive_power_limits(device)
        if (q_set_point >= q_limits.max + BOUNDS_TOLERANCE) ||
           (q_set_point <= q_limits.min - BOUNDS_TOLERANCE)
            @error "Unit $(PSY.get_name(device)) Q=$(q_set_point) above limits. Q_max = $(q_limits.max) Q_min = $(q_limits.min)"
        end
    end

    return
end

"""
Updates system voltages and powers with power flow results
"""
function write_powerflow_solution!(sys::PSY.System, result::Vector{Float64})
    buses = enumerate(
        sort!(collect(PSY.get_components(PSY.Bus, sys)); by = x -> PSY.get_number(x)),
    )

    for (ix, bus) in buses
        if bus.bustype == PSY.BusTypes.REF
            P_gen = result[2 * ix - 1]
            Q_gen = result[2 * ix]
            _power_redistribution_ref(sys, P_gen, Q_gen, bus)
        elseif bus.bustype == PSY.BusTypes.PV
            Q_gen = result[2 * ix - 1]
            bus.angle = result[2 * ix]
            _reactive_power_redistribution_pv(sys, Q_gen, bus)
        elseif bus.bustype == PSY.BusTypes.PQ
            Vm = result[2 * ix - 1]
            θ = result[2 * ix]
            PSY.set_magnitude!(bus, Vm)
            PSY.set_angle!(bus, θ)
        end
    end

    _update_branch_flow!(sys)
    return
end

# returns list of branches names and buses numbers: ABA case
function _get_branches_buses(data::ABAPowerFlowData)
    return axes(data.aux_network_matrix)[1], axes(data.aux_network_matrix)[2]
end

# returns list of branches names and buses numbers: PTDF case
function _get_branches_buses(data::PTDFPowerFlowData)
    return axes(data.power_network_matrix)[2], axes(data.power_network_matrix)[1]
end

# ! missing version for virtual PTDF

"""
Return power flow results in dictionary of dataframes.
"""
function write_results(
    data::Union{PTDFPowerFlowData, ABAPowerFlowData},
    sys::PSY.System
    )

    @info("Voltages are exported in pu. Powers are exported in MW/MVAr.")
    
    # get bus and branches
    branches, buses = _get_branches_buses(data)

    # get branches from/to buses
    from_bus = Vector{Int}(undef, length(branches))
    to_bus = Vector{Int}(undef, length(branches))
    for (i, branch) in enumerate(branches)
        br = PSY.get_component(PSY.ACBranch, sys, branch)
        from_bus[i] = PSY.get_number(PSY.get_arc(br).from)
        to_bus[i] = PSY.get_number(PSY.get_arc(br).to)
    end

    # get injections and withdrawals
    P_gen_vect = data.bus_activepower_injection
    Q_gen_vect = data.bus_reactivepower_injection
    P_load_vect = data.bus_activepower_withdrawals
    Q_load_vect = data.bus_reactivepower_withdrawals

    # get flows on each line
    P_from_to_vect = zeros(length(branches))
    Q_from_to_vect = zeros(length(branches))
    P_to_from_vect = zeros(length(branches))
    Q_to_from_vect = zeros(length(branches))
    for i in 1:length(branches)
        if data.branch_flow_values[i] >= 0
            P_from_to_vect[i] = data.branch_flow_values[i]
            P_to_from_vect[i] = 0
        else
            P_from_to_vect[i] = 0
            P_to_from_vect[i] = data.branch_flow_values[i]
        end
    end

    bus_df = DataFrames.DataFrame(;
        bus_number = buses,
        Vm = data.bus_magnitude,
        θ = data.bus_angles,
        P_gen = P_gen_vect,
        P_load = P_load_vect,
        P_net = P_gen_vect - P_load_vect,
        Q_gen = Q_gen_vect,
        Q_load = Q_load_vect,
        Q_net = Q_gen_vect - Q_load_vect,
    )
    DataFrames.sort!(bus_df, :bus_number)

    branch_df = DataFrames.DataFrame(;
        line_name = branches,
        bus_from = from_bus,
        bus_to = to_bus,
        P_from_to = P_from_to_vect,
        Q_from_to = Q_from_to_vect,
        P_to_from = P_to_from_vect,
        Q_to_from = Q_to_from_vect,
        P_losses = zeros(length(branches)),
        Q_losses = zeros(length(branches)),
    )
    DataFrames.sort!(branch_df, [:bus_from, :bus_to])

    return Dict("bus_results" => bus_df, "flow_results" => branch_df)
end

function write_results(
    ::ACPowerFlow,
    sys::PSY.System,
    result::Vector{Float64}
    )

    @info("Voltages are exported in pu. Powers are exported in MW/MVAr.")
    buses = sort!(collect(PSY.get_components(PSY.Bus, sys)); by = x -> PSY.get_number(x))
    N_BUS = length(buses)
    bus_map = Dict(buses .=> 1:N_BUS)
    sys_basepower = PSY.get_base_power(sys)
    sources = PSY.get_components(d -> !isa(d, PSY.ElectricLoad), PSY.StaticInjection, sys)
    Vm_vect = fill(0.0, N_BUS)
    θ_vect = fill(0.0, N_BUS)
    P_gen_vect = fill(0.0, N_BUS)
    Q_gen_vect = fill(0.0, N_BUS)
    P_load_vect = fill(0.0, N_BUS)
    Q_load_vect = fill(0.0, N_BUS)

    for (ix, bus) in enumerate(buses)
        P_load_vect[ix], Q_load_vect[ix] = _get_load_data(sys, bus) .* sys_basepower
        P_admittance, Q_admittance = _get_fixed_admittance_power(sys, bus, result, ix)
        P_load_vect[ix] += P_admittance * sys_basepower
        Q_load_vect[ix] += Q_admittance * sys_basepower
        if bus.bustype == PSY.BusTypes.REF
            Vm_vect[ix] = PSY.get_magnitude(bus)
            θ_vect[ix] = PSY.get_angle(bus)
            P_gen_vect[ix] = result[2 * ix - 1] * sys_basepower
            Q_gen_vect[ix] = result[2 * ix] * sys_basepower
        elseif bus.bustype == PSY.BusTypes.PV
            Vm_vect[ix] = PSY.get_magnitude(bus)
            θ_vect[ix] = result[2 * ix]
            for gen in sources
                !PSY.get_available(gen) && continue
                if gen.bus == bus
                    P_gen_vect[ix] += PSY.get_active_power(gen) * sys_basepower
                end
            end
            Q_gen_vect[ix] = result[2 * ix - 1] * sys_basepower
        elseif bus.bustype == PSY.BusTypes.PQ
            Vm_vect[ix] = result[2 * ix - 1]
            θ_vect[ix] = result[2 * ix]
            for gen in sources
                !PSY.get_available(gen) && continue
                if gen.bus == bus
                    P_gen_vect[ix] += PSY.get_active_power(gen) * sys_basepower
                    Q_gen_vect[ix] += PSY.get_reactive_power(gen) * sys_basepower
                end
            end
        end
    end

    branches = PSY.get_components(PSY.ACBranch, sys)
    N_BRANCH = length(branches)
    P_from_to_vect = fill(0.0, N_BRANCH)
    Q_from_to_vect = fill(0.0, N_BRANCH)
    P_to_from_vect = fill(0.0, N_BRANCH)
    Q_to_from_vect = fill(0.0, N_BRANCH)
    for (ix, b) in enumerate(branches)
        !PSY.get_available(b) && continue
        bus_f_ix = bus_map[PSY.get_arc(b).from]
        bus_t_ix = bus_map[PSY.get_arc(b).to]
        V_from = Vm_vect[bus_f_ix] * (cos(θ_vect[bus_f_ix]) + sin(θ_vect[bus_f_ix]) * 1im)
        V_to = Vm_vect[bus_t_ix] * (cos(θ_vect[bus_t_ix]) + sin(θ_vect[bus_t_ix]) * 1im)
        P_from_to_vect[ix], Q_from_to_vect[ix] = flow_func(b, V_from, V_to) .* sys_basepower
        P_to_from_vect[ix], Q_to_from_vect[ix] = flow_func(b, V_to, V_from) .* sys_basepower
    end

    bus_df = DataFrames.DataFrame(;
        bus_number = PSY.get_number.(buses),
        Vm = Vm_vect,
        θ = θ_vect,
        P_gen = P_gen_vect,
        P_load = P_load_vect,
        P_net = P_gen_vect - P_load_vect,
        Q_gen = Q_gen_vect,
        Q_load = Q_load_vect,
        Q_net = Q_gen_vect - Q_load_vect,
    )

    branch_df = DataFrames.DataFrame(;
        line_name = PSY.get_name.(branches),
        bus_from = PSY.get_number.(PSY.get_from.(PSY.get_arc.(branches))),
        bus_to = PSY.get_number.(PSY.get_to.(PSY.get_arc.(branches))),
        P_from_to = P_from_to_vect,
        Q_from_to = Q_from_to_vect,
        P_to_from = P_to_from_vect,
        Q_to_from = Q_to_from_vect,
        P_losses = P_from_to_vect + P_to_from_vect,
        Q_losses = Q_from_to_vect + Q_to_from_vect,
    )
    DataFrames.sort!(branch_df, [:bus_from, :bus_to])

    return Dict("bus_results" => bus_df, "flow_results" => branch_df)
end