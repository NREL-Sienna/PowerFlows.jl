_SingleComponentLoad = Union{PSY.PowerLoad, PSY.ExponentialLoad, PSY.InterruptiblePowerLoad}

"""
Return the reactive power limits that should be used in power flow calculations and PSS/E
exports. Redirects to `PSY.get_reactive_power_limits` in all but special cases.
"""
get_reactive_power_limits_for_power_flow(gen::PSY.Device) =
    PSY.get_reactive_power_limits(gen)

function get_reactive_power_limits_for_power_flow(gen::PSY.RenewableNonDispatch)
    val = PSY.get_reactive_power(gen)
    return (min = val, max = val)
end

function get_reactive_power_limits_for_power_flow(gen::PSY.Storage)
    limits = PSY.get_reactive_power_limits(gen)
    isnothing(limits) && return (min = -Inf, max = Inf)  # TODO decide on proper behavior in this case
    return limits
end

"""
Return the active power limits that should be used in power flow calculations and PSS/E
exports. Redirects to `PSY.get_active_power_limits` in all but special cases.
"""
get_active_power_limits_for_power_flow(gen::PSY.Device) = PSY.get_active_power_limits(gen)

get_active_power_limits_for_power_flow(::PSY.Source) = (min = -Inf, max = Inf)

function get_active_power_limits_for_power_flow(gen::PSY.RenewableNonDispatch)
    val = PSY.get_active_power(gen)
    return (min = val, max = val)
end

get_active_power_limits_for_power_flow(gen::PSY.RenewableDispatch) =
    (min = 0.0, max = PSY.get_rating(gen))

# TODO verify whether this is the correct behavior for Storage, (a) for redistribution and (b) for exporting
get_active_power_limits_for_power_flow(gen::PSY.Storage) =
    (min = 0.0, max = PSY.get_output_active_power_limits(gen).max)

_get_bus_ix(
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    bus_number::Int,
) =
    bus_lookup[get(reverse_bus_search_map, bus_number, bus_number)]

function _get_injections!(
    bus_activepower_injection::Vector{Float64},
    bus_reactivepower_injection::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    sys::PSY.System,
)
    for source in PSY.get_available_components(PSY.StaticInjection, sys)
        if contributes_active_power(source) &&
           active_power_contribution_type(source) == PowerContributionType.INJECTION
            bus = PSY.get_bus(source)
            bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
            bus_activepower_injection[bus_ix] += PSY.get_active_power(source)
        end
        if contributes_reactive_power(source) &&
           reactive_power_contribution_type(source) == PowerContributionType.INJECTION
            bus = PSY.get_bus(source)
            bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
            # yet to implement control mode etc. for FACTS devices.
            if source isa PSY.FACTSControlDevice
                bus_reactivepower_injection[bus_ix] +=
                    PSY.get_reactive_power_required(source)
            else
                bus_reactivepower_injection[bus_ix] += PSY.get_reactive_power(source)
            end
        end
    end
    return
end

function _get_withdrawals!(
    bus_activepower_withdrawals::Vector{Float64},
    bus_reactivepower_withdrawals::Vector{Float64},
    bus_activepower_constant_current_withdrawals::Vector{Float64},
    bus_reactivepower_constant_current_withdrawals::Vector{Float64},
    bus_activepower_constant_impedance_withdrawals::Vector{Float64},
    bus_reactivepower_constant_impedance_withdrawals::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    sys::PSY.System,
)
    # constant power withdrawals
    for l in PSY.get_available_components(PSY.StaticInjection, sys)
        if contributes_active_power(l) &&
           active_power_contribution_type(l) == PowerContributionType.WITHDRAWAL
            bus = PSY.get_bus(l)
            bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
            bus_activepower_withdrawals[bus_ix] += PSY.get_active_power(l)
        end
        if contributes_reactive_power(l) &&
           reactive_power_contribution_type(l) == PowerContributionType.WITHDRAWAL
            bus = PSY.get_bus(l)
            bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
            bus_reactivepower_withdrawals[bus_ix] += PSY.get_reactive_power(l)
        end
    end
    # handle StandardLoad: they have constant current and constant impedance withdrawals,
    # and the getter for constant power is named differently (get_constant_active_power)
    for l in PSY.get_available_components(
        Union{PSY.StandardLoad, PSY.InterruptibleStandardLoad},
        sys,
    )
        bus = PSY.get_bus(l)
        bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
        bus_activepower_withdrawals[bus_ix] += PSY.get_constant_active_power(l)
        bus_reactivepower_withdrawals[bus_ix] += PSY.get_constant_reactive_power(l)
        bus_activepower_constant_current_withdrawals[bus_ix] +=
            PSY.get_current_active_power(l)
        bus_activepower_constant_impedance_withdrawals[bus_ix] +=
            PSY.get_impedance_active_power(l)
        bus_reactivepower_constant_current_withdrawals[bus_ix] +=
            PSY.get_current_reactive_power(l)
        bus_reactivepower_constant_impedance_withdrawals[bus_ix] +=
            PSY.get_impedance_reactive_power(l)
    end
    # FixedAdmittance components are already included in the Ybus matrix.
    for sa in PSY.get_available_components(PSY.SwitchedAdmittance, sys)
        bus = PSY.get_bus(sa)
        bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
        Y = PSY.get_Y(sa) + sum(PSY.get_initial_status(sa) .* PSY.get_Y_increase(sa))
        # Here we implement the switched admittance element as a constant impedance load.
        # The inputs for ZIP loads are provided for V = 1.0 p.u., so
        # the following is equivalent to S = V * conj(Y * V) for V = 1.0 p.u.
        # As conj(Y) = G - jB, we have P = G and Q = -B
        # (we could use +=real(conj(Y)) and +=imag(conj(Y)) as well).
        bus_activepower_constant_impedance_withdrawals[bus_ix] += real(Y)
        bus_reactivepower_constant_impedance_withdrawals[bus_ix] -= imag(Y)
    end
    for sc in PSY.get_available_components(PSY.SynchronousCondenser, sys)
        bus = PSY.get_bus(sc)
        bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
        bus_activepower_withdrawals[bus_ix] += PSY.get_active_power_losses(sc)
        # reactive power handled already:
        # contributes_reactive_power(PSY.SynchronousCondenser) is true.
    end
    return
end

# TODO: Might need changes if we have SwitchedAdmittances
function _get_reactive_power_bound!(
    bus_reactivepower_bounds::Vector{Tuple{Float64, Float64}},
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    sys::PSY.System,
)
    for source in PSY.get_available_components(PSY.StaticInjection, sys)
        isa(source, PSY.ElectricLoad) && continue
        isa(source, PSY.FACTSControlDevice) && continue # FIXME: FACTS devices.
        bus = PSY.get_bus(source)
        bus_ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
        reactive_power_limits = get_reactive_power_limits_for_power_flow(source)
        if reactive_power_limits !== nothing
            bus_reactivepower_bounds[bus_ix] = (
                bus_reactivepower_bounds[bus_ix][1] + min(0, reactive_power_limits.min),
                bus_reactivepower_bounds[bus_ix][2] + max(0, reactive_power_limits.max),
            )
        else
            @warn("Reactive Power Bounds at Bus $(PSY.get_name(bus)) set to (-Inf, Inf)")
            bus_reactivepower_bounds[bus_ix] = (-Inf, Inf)
        end
    end
    return
end

function _initialize_bus_data!(
    bus_type::Vector{PSY.ACBusTypes},
    bus_angles::Vector{Float64},
    bus_magnitude::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    bus_reduction_map::Dict{Int, Set{Int}},
    reverse_bus_search_map::Dict{Int, Int},
    sys::PSY.System,
    correct_bustypes::Bool = false,
)
    # correct/validate the bus types.
    forced_PV = must_be_PV(sys)
    possible_PV = can_be_PV(sys)
    bus_numbers = PSY.get_bus_numbers(sys)
    temp_bus_types = Dict{Int, PSY.ACBusTypes}()
    sizehint!(temp_bus_types, length(bus_numbers))
    temp_bus_map = Dict{Int, String}()
    sizehint!(temp_bus_map, length(bus_numbers))
    for bus in PSY.get_components(PSY.ACBus, sys)
        bt = PSY.get_bustype(bus)
        bus_no = PSY.get_number(bus)
        bus_name = PSY.get_name(bus)
        temp_bus_map[bus_no] = bus_name
        if (bt == PSY.ACBusTypes.PV || bt == PSY.ACBusTypes.REF) && !(bus_no in possible_PV)
            if correct_bustypes
                @warn "No available sources at bus $bus_name  of bus type 2 (PV). " *
                      "Treating that bus as PQ for purposes of the power flow." maxlog =
                    PF_MAX_LOG
                bt = PSY.ACBusTypes.PQ
            else
                throw(
                    ArgumentError(
                        "No available sources at bus $bus_name of bus type 2 (PV)." *
                        " Please change the bus type to PQ.",
                    ),
                )
            end
        elseif bt == PSY.ACBusTypes.PQ && bus_no in forced_PV
            @warn "Active generators found at bus $bus_name of bus type 1 (PQ), i.e. " *
                  "different than 2 (PV). Consider checking your data inputs." maxlog =
                PF_MAX_LOG
        end
        temp_bus_types[bus_no] = bt
    end

    # FIXME handle combining buses more intelligently.
    # perhaps move some of this logic into PNM
    for (bus_no, reduced_bus_nos) in bus_reduction_map
        # pick the "highest" bus type among the reduced buses, where REF > PV > PQ.
        corrected_bus_types =
            [temp_bus_types[reduced_bus_no] for reduced_bus_no in reduced_bus_nos]
        push!(corrected_bus_types, temp_bus_types[bus_no])
        combined_bus_type = findmax(bt -> BUS_TYPE_PRIORITIES[bt], corrected_bus_types)[1]
        ix = _get_bus_ix(bus_lookup, reverse_bus_search_map, bus_no)
        bus_type[ix] = combined_bus_type
        # TODO: combine the angles/magnitudes, between the reduced buses?
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.ACBus, sys, bus_name)
        bus_angles[ix] = PSY.get_angle(bus)
        bus_vm = PSY.get_magnitude(bus)
        # prevent unfeasible starting values for voltage magnitude at PQ buses (for PV and REF buses we cannot do this):
        if combined_bus_type == PSY.ACBusTypes.PQ &&
           bus_vm < BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN
            @warn(
                "Initial bus voltage magnitude of $bus_vm p.u. at PQ bus $bus_name is below the plausible minimum cut-off value of $BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN p.u. and has been set to $BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN p.u.",
                maxlog = PF_MAX_LOG,
            )
            bus_vm = BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN
        elseif combined_bus_type == PSY.ACBusTypes.PQ &&
               bus_vm > BUS_VOLTAGE_MAGNITUDE_CUTOFF_MAX
            @warn(
                "Initial bus voltage magnitude of $bus_vm p.u. at PQ bus $bus_name is above the plausible maximum cut-off value of $BUS_VOLTAGE_MAGNITUDE_CUTOFF_MAX p.u. and has been set to $BUS_VOLTAGE_MAGNITUDE_CUTOFF_MAX p.u.",
                maxlog = PF_MAX_LOG,
            )
            bus_vm = BUS_VOLTAGE_MAGNITUDE_CUTOFF_MAX
        end
        bus_magnitude[ix] = bus_vm
    end
end
##############################################################################
# Matrix Methods #############################################################

"""Matrix multiplication A*x. Written this way because a VirtualPTDF 
matrix does not store all of its entries: instead, it calculates
them (or retrieves them from cache), one element or one row at a time."""
function my_mul_mt(
    A::PNM.VirtualPTDF,
    x::Vector{Float64},
)
    y = zeros(length(A.axes[1]))
    for i in 1:length(A.axes[1])
        name_ = A.axes[1][i]
        y[i] = dot(A[name_, :], x)
    end
    return y
end

"""Similar to above: A*X where X is a matrix."""
my_mul_mt(
    A::PNM.VirtualPTDF,
    X::Matrix{Float64},
) = vcat((A[name_, :]' * X for name_ in A.axes[1])...)

function make_dc_powerflowdata(
    sys,
    time_steps,
    timestep_names,
    power_network_matrix,
    aux_network_matrix,
    bus_lookup,
    arc_lookup,
    valid_ix,
    converged,
    loss_factors,
    calculate_loss_factors,
    voltage_stability_factors,
    calculate_voltage_stability_factors,
    generator_slack_participation_factors,
    correct_bustypes,
)
    timestep_map = Dict(zip([i for i in 1:time_steps], timestep_names))
    neighbors = Vector{Set{Int}}()
    return make_powerflowdata(
        sys,
        time_steps,
        power_network_matrix,
        aux_network_matrix,
        bus_lookup,
        arc_lookup,
        timestep_map,
        valid_ix,
        neighbors,
        converged,
        loss_factors,
        calculate_loss_factors,
        voltage_stability_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        correct_bustypes,
    )
end

function _add_gspf_to_ijv!(
    I::Vector{Int},
    J::Vector{Int},
    V::Vector{Float64},
    sys::System,
    gsp_factors::Dict{Tuple{DataType, String}, Float64},
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    time_steps_iter::AbstractVector{<:Integer},
)
    for ((gen_type, gen_name), val) in gsp_factors
        val == 0.0 && continue
        gen = PSY.get_component(gen_type, sys, gen_name)
        isnothing(gen) && throw(ArgumentError("$gen_type $gen_name not found"))
        PSY.get_available(gen) || continue
        bus = PSY.get_bus(gen)
        bus_idx = _get_bus_ix(bus_lookup, reverse_bus_search_map, PSY.get_number(bus))
        for time_step in time_steps_iter
            push!(I, bus_idx)
            push!(J, time_step)
            push!(V, val)
        end
    end
    return
end

function make_bus_slack_participation_factors(
    sys::System,
    generator_slack_participation_factors::Dict{Tuple{DataType, String}, Float64},
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    time_steps::Int,
    n_buses::Int,
    ::Matrix{PSY.ACBusTypes},
)
    I = Int[]
    J = Int[]
    V = Float64[]

    _add_gspf_to_ijv!(
        I,
        J,
        V,
        sys,
        generator_slack_participation_factors,
        bus_lookup,
        reverse_bus_search_map,
        1:time_steps,
    )

    bus_slack_participation_factors = sparse(I, J, V, n_buses, time_steps)
    return bus_slack_participation_factors,
    repeat([generator_slack_participation_factors], time_steps)
end

function make_bus_slack_participation_factors(
    sys::System,
    generator_slack_participation_factors::Vector{Dict{Tuple{DataType, String}, Float64}},
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    time_steps::Int,
    n_buses::Int,
    bus_type::Matrix{PSY.ACBusTypes},
)
    if length(generator_slack_participation_factors) == 1
        return make_bus_slack_participation_factors(
            sys,
            generator_slack_participation_factors[1],
            bus_lookup,
            reverse_bus_search_map,
            time_steps,
            n_buses,
            bus_type,
        )
    end

    if length(generator_slack_participation_factors) < time_steps
        throw(
            ArgumentError(
                "slack_participation_factors must have at least the same length as time_steps",
            ),
        )
    end

    # A sparse matrix constructor is used here, and the duplicates at the same locations are summed by default.
    # This way, the generator slack participation factors are aggregated per bus.
    I = Int[]
    J = Int[]
    V = Float64[]

    for (time_step, factors) in enumerate(generator_slack_participation_factors)
        _add_gspf_to_ijv!(
            I,
            J,
            V,
            sys,
            factors,
            bus_lookup,
            reverse_bus_search_map,
            time_step:time_step,
        )
    end

    bus_slack_participation_factors = sparse(I, J, V, n_buses, time_steps)
    return bus_slack_participation_factors, generator_slack_participation_factors
end

function make_bus_slack_participation_factors(
    ::System,
    ::Nothing,
    ::Dict{Int, Int},
    ::Dict{Int, Int},
    time_steps::Int,
    n_buses::Int,
    bus_type::Matrix{PSY.ACBusTypes},
)
    I = Int[]
    J = Int[]
    V = Float64[]

    for time_step in 1:time_steps
        for (ix, bt) in enumerate(bus_type[:, time_step])
            bt == PSY.ACBusTypes.REF || continue
            push!(I, ix)
            push!(J, time_step)
            push!(V, 1.0)
        end
    end

    bus_slack_participation_factors = sparse(I, J, V, n_buses, time_steps)

    return bus_slack_participation_factors, nothing
end

"""Return set of all bus numbers that must be PV: i.e. have an available generator."""
function must_be_PV(sys::System)
    gen_buses = Set{Int}()
    for gen in PSY.get_available_components(PSY.Generator, sys)
        push!(gen_buses, PSY.get_number(PSY.get_bus(gen)))
    end
    # PSSe counts buses with switched shunts as PV, so we do the same here.
    for gen in PSY.get_available_components(PSY.SwitchedAdmittance, sys)
        push!(gen_buses, PSY.get_number(PSY.get_bus(gen)))
    end
    for gen in PSY.get_available_components(PSY.SynchronousCondenser, sys)
        push!(gen_buses, PSY.get_number(PSY.get_bus(gen)))
    end
    return gen_buses
end

"""Return set of all bus numbers that can be PV: i.e. have an available generator,
or certain voltage regulation devices."""
function can_be_PV(sys::System)
    source_buses = must_be_PV(sys)
    for source in PSY.get_available_components(PSY.Source, sys)
        push!(source_buses, PSY.get_number(PSY.get_bus(source)))
    end
    return source_buses
end

function make_powerflowdata(
    sys,
    time_steps,
    power_network_matrix,
    aux_network_matrix,
    bus_lookup,
    arc_lookup,
    timestep_map,
    valid_ix,
    neighbors,
    converged,
    loss_factors,
    calculate_loss_factors,
    voltage_stability_factors = nothing,
    calculate_voltage_stability_factors = nothing,
    generator_slack_participation_factors = nothing,
    correct_bustypes::Bool = false,
    lcc_p_set = nothing,
)
    n_buses = length(bus_lookup)
    n_arcs = length(arc_lookup)
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = ones(Float64, n_buses)
    nrd = PNM.get_network_reduction_data(power_network_matrix)
    reverse_bus_search_map = PNM.get_reverse_bus_search_map(nrd)
    bus_reduction_map = PNM.get_bus_reduction_map(nrd)

    _initialize_bus_data!(
        bus_type,
        bus_angles,
        bus_magnitude,
        bus_lookup,
        bus_reduction_map,
        reverse_bus_search_map,
        sys,
        correct_bustypes,
    )

    # define injection vectors related to the first timestep
    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    _get_injections!(
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_lookup,
        reverse_bus_search_map,
        sys,
    )

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    bus_activepower_constant_current_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_constant_current_withdrawals = zeros(Float64, n_buses)
    bus_activepower_constant_impedance_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_constant_impedance_withdrawals = zeros(Float64, n_buses)
    _get_withdrawals!(
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_activepower_constant_current_withdrawals,
        bus_reactivepower_constant_current_withdrawals,
        bus_activepower_constant_impedance_withdrawals,
        bus_reactivepower_constant_impedance_withdrawals,
        bus_lookup,
        reverse_bus_search_map,
        sys,
    )

    # Define fields as matrices whose number of columns is equal to the number of time_steps
    bus_activepower_injection_1 = zeros(Float64, n_buses, time_steps)
    bus_reactivepower_injection_1 = zeros(Float64, n_buses, time_steps)
    bus_activepower_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    bus_reactivepower_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    bus_activepower_constant_current_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    bus_reactivepower_constant_current_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    bus_activepower_constant_impedance_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    bus_reactivepower_constant_impedance_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    # TODO: our bus reactive power bounds are undefined for all timesteps besides the first?!
    bus_reactivepower_bounds_1 = fill((NaN, NaN), n_buses, time_steps)
    bus_magnitude_1 = ones(Float64, n_buses, time_steps)
    bus_angles_1 = zeros(Float64, n_buses, time_steps)

    # Initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_activepower_constant_current_withdrawals_1[:, 1] .=
        bus_activepower_constant_current_withdrawals
    bus_reactivepower_constant_current_withdrawals_1[:, 1] .=
        bus_reactivepower_constant_current_withdrawals
    bus_activepower_constant_impedance_withdrawals_1[:, 1] .=
        bus_activepower_constant_impedance_withdrawals
    bus_reactivepower_constant_impedance_withdrawals_1[:, 1] .=
        bus_reactivepower_constant_impedance_withdrawals
    bus_magnitude_1[:, 1] .= bus_magnitude
    bus_angles_1[:, 1] .= bus_angles

    bus_reactivepower_bounds = Vector{Tuple{Float64, Float64}}(undef, n_buses)
    for i in 1:n_buses
        bus_reactivepower_bounds[i] = (0.0, 0.0)
    end
    _get_reactive_power_bound!(
        bus_reactivepower_bounds,
        bus_lookup,
        reverse_bus_search_map,
        sys,
    )
    bus_reactivepower_bounds_1[:, 1] .= bus_reactivepower_bounds

    # Initial bus types are same for every time period
    bus_type_1 = repeat(bus_type; outer = [1, time_steps])
    @assert size(bus_type_1) == (n_buses, time_steps)

    # Initial slack participation factors are same for every time period
    bus_slack_participation_factors, generator_slack_participation_factors =
        make_bus_slack_participation_factors(
            sys,
            generator_slack_participation_factors,
            bus_lookup,
            reverse_bus_search_map,
            time_steps,
            n_buses,
            bus_type_1,
        )

    # Initial flows are all zero
    arc_activepower_flow_from_to = zeros(Float64, n_arcs, time_steps)
    arc_reactivepower_flow_from_to = zeros(Float64, n_arcs, time_steps)
    arc_activepower_flow_to_from = zeros(Float64, n_arcs, time_steps)
    arc_reactivepower_flow_to_from = zeros(Float64, n_arcs, time_steps)

    lcc = LCCParameters(sys, bus_lookup, reverse_bus_search_map, time_steps)

    return PowerFlowData(
        bus_lookup,
        arc_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        bus_activepower_constant_current_withdrawals_1,
        bus_reactivepower_constant_current_withdrawals_1,
        bus_activepower_constant_impedance_withdrawals_1,
        bus_reactivepower_constant_impedance_withdrawals_1,
        bus_reactivepower_bounds_1,
        generator_slack_participation_factors,
        bus_slack_participation_factors,
        bus_type_1,
        bus_magnitude_1,
        bus_angles_1,
        arc_activepower_flow_from_to,
        arc_reactivepower_flow_from_to,
        arc_activepower_flow_to_from,
        arc_reactivepower_flow_to_from,
        timestep_map,
        valid_ix,
        power_network_matrix,
        aux_network_matrix,
        neighbors,
        converged,
        loss_factors,
        calculate_loss_factors,
        voltage_stability_factors,
        calculate_voltage_stability_factors,
        lcc,
    )
end

function validate_voltages(x::Vector{Float64},
    bus_types::AbstractArray{PSY.ACBusTypes},
    range::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = VM_VALIDATION_RANGE,
    i::Int64 = 1,
)
    # outside_range = sizehint!(Vector{Int64}(), MAX_INDS_TO_PRINT)
    num_outside_range = 0
    for (i, bt) in enumerate(bus_types)
        if bt == PSY.ACBusTypes.PQ
            if (x[2 * i - 1] < range.min || x[2 * i - 1] > range.max) #&&
                # size(outside_range, 1) <= MAX_INDS_TO_PRINT
                #push!(outside_range, i)
                num_outside_range += 1
            end
        end
    end
    if num_outside_range > 0
        @warn "Iteration $i: voltage magnitudes outside of range $range at " *
              "$num_outside_range PQ buses." maxlog = PF_MAX_LOG
    end
    #=
    if size(outside_range, 1) > MAX_INDS_TO_PRINT
        @warn "Iteration $i: voltage magnitudes outside of range $range at over $MAX_INDS_TO_PRINT buses." maxlog =
            PF_MAX_LOG
    elseif size(outside_range, 1) > 0
        @warn "Iteration $i: voltage magnitudes outside of range $range at $(size(outside_range, 1)) buses: $(outside_range)" maxlog =
            PF_MAX_LOG
    end
    =#
    return
end

"""Weighted dot product of two vectors."""
wdot(wx::Vector{Float64}, x::Vector{Float64}, wy::Vector{Float64}, y::Vector{Float64}) =
    LinearAlgebra.dot(wx .* x, wy .* y)

"""Weighted norm of two vectors."""
wnorm(w::Vector{Float64}, x::Vector{Float64}) = norm(w .* x)
"""For pretty printing floats in debugging messages."""
siground(x::Float64) = round(x; sigdigits = 3)

"""
    _calculate_ϕ_lcc(α::Float64, I_dc::Float64, x_t::Float64, Vm::Float64) -> Float64

Compute the phase angle ϕ for LCC converter calculations.
"""
function _calculate_ϕ_lcc(
    t::Float64,
    α::Float64,
    I_dc::Float64,
    x_t::Float64,
    Vm::Float64,
)::Float64
    return acos(clamp(sign(I_dc) * (cos(α) - (x_t * I_dc) / (sqrt(2) * Vm * t)), -1.0, 1.0))
end

"""
    _calculate_y_lcc(t::Float64, I_dc::Float64, Vm::Float64, ϕ::Float64) -> ComplexF64

Compute the admittance value Y for LCC converter calculations.
"""
function _calculate_y_lcc(t::Float64, I_dc::Float64, Vm::Float64, ϕ::Float64)::ComplexF64
    return t / Vm * sqrt(6) / π * I_dc * exp(-1im * ϕ)
end

"""
    _calculate_dQ_dV_lcc(t::Float64, I_dc::Float64, x_t::Float64, Vm::Float64, ϕ::Float64) -> Float64

Compute the derivative of reactive power Q with respect to voltage magnitude Vm for LCC converter calculations.
"""
function _calculate_dQ_dV_lcc(
    t::Float64,
    I_dc::Float64,
    x_t::Float64,
    Vm::Float64,
    ϕ::Float64,
)::Float64
    return Vm * t * sqrt(6) / π * I_dc *
           (sin(ϕ) - cos(ϕ) * x_t * I_dc / (sqrt(2) * Vm * t * sin(ϕ)^2))
end

"""
    _calculate_dQ_dt_lcc(t::Float64, I_dc::Float64, x_t::Float64, Vm::Float64, ϕ::Float64) -> Float64

Compute the derivative of reactive power Q with respect to transformer tap t for LCC converter calculations.
"""
function _calculate_dQ_dt_lcc(
    t::Float64,
    I_dc::Float64,
    x_t::Float64,
    Vm::Float64,
    ϕ::Float64,
)::Float64
    return Vm * t * sqrt(6) / π * I_dc *
           (sin(ϕ) / t - cos(ϕ) * x_t * I_dc / (sqrt(2) * Vm * t^2 * sin(ϕ)^2))
end

"""
    _calculate_dQ_dα_lcc(t::Float64, I_dc::Float64, x_t::Float64, Vm::Float64, ϕ::Float64, α::Float64) -> Float64

Compute the derivative of reactive power Q with respect to firing/extinction angle α for LCC converter calculations.
"""
function _calculate_dQ_dα_lcc(
    t::Float64,
    I_dc::Float64,
    x_t::Float64,
    Vm::Float64,
    ϕ::Float64,
    α::Float64,
)::Float64
    return Vm * t * sqrt(6) / π * I_dc * cos(ϕ) * sin(α) / sqrt(1 - cos(ϕ)^2)
end

function _update_ybus_lcc!(ybus_facts, data, time_step)
    for i in eachindex(data.lcc.rectifier.bus)
        data.lcc.rectifier.phi[i, time_step] = _calculate_ϕ_lcc(
            data.lcc.rectifier.tap[i, time_step],
            data.lcc.rectifier.thyristor_angle[i, time_step],
            data.lcc.i_dc[i, time_step],
            data.lcc.rectifier.transformer_reactance[i],
            data.bus_magnitude[data.lcc.rectifier.bus[i], time_step],
        )
        data.lcc.inverter.phi[i, time_step] = _calculate_ϕ_lcc(
            data.lcc.inverter.tap[i, time_step],
            data.lcc.inverter.thyristor_angle[i, time_step],
            -data.lcc.i_dc[i, time_step],
            data.lcc.inverter.transformer_reactance[i],
            data.bus_magnitude[data.lcc.inverter.bus[i], time_step],
        )
        ybus_facts[data.lcc.rectifier.bus[i], data.lcc.rectifier.bus[i]] = _calculate_y_lcc(
            data.lcc.rectifier.tap[i, time_step],
            data.lcc.i_dc[i, time_step],
            data.bus_magnitude[data.lcc.rectifier.bus[i], time_step],
            data.lcc.rectifier.phi[i, time_step],
        )
        ybus_facts[data.lcc.inverter.bus[i], data.lcc.inverter.bus[i]] = _calculate_y_lcc(
            data.lcc.inverter.tap[i, time_step],
            data.lcc.i_dc[i, time_step],
            data.bus_magnitude[data.lcc.inverter.bus[i], time_step],
            data.lcc.inverter.phi[i, time_step],
        )
        data.power_network_matrix.branch_admittance_from_to[
            data.lcc.arc[i],
            data.lcc.arc[i].from,
        ] = ybus_facts[data.lcc.rectifier.bus[i], data.lcc.rectifier.bus[i]]
        data.power_network_matrix.branch_admittance_to_from[
            data.lcc.arc[i],
            data.lcc.arc[i].to,
        ] = ybus_facts[data.lcc.inverter.bus[i], data.lcc.inverter.bus[i]]
    end
    return
end
