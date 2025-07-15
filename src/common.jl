_SingleComponentLoad = Union{PSY.PowerLoad, PSY.ExponentialLoad, PSY.InterruptiblePowerLoad}
get_total_p(l::_SingleComponentLoad) = PSY.get_active_power(l)
get_total_q(l::_SingleComponentLoad) = PSY.get_reactive_power(l)

function get_total_p(l::PSY.StandardLoad)
    return PSY.get_constant_active_power(l) +
           PSY.get_current_active_power(l) +
           PSY.get_impedance_active_power(l)
end

function get_total_q(l::PSY.StandardLoad)
    return PSY.get_constant_reactive_power(l) +
           PSY.get_current_reactive_power(l) +
           PSY.get_impedance_reactive_power(l)
end

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

function _get_injections!(
    bus_activepower_injection::Vector{Float64},
    bus_reactivepower_injection::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System,
)
    for source in PSY.get_components(PSY.StaticInjection, sys)
        isa(source, PSY.ElectricLoad) && continue
        !PSY.get_available(source) && continue
        bus = PSY.get_bus(source)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        # see issue #1463 in PSY
        if !isa(source, PSY.SynchronousCondenser)
            bus_activepower_injection[bus_ix] += PSY.get_active_power(source)
        end
        bus_reactivepower_injection[bus_ix] += PSY.get_reactive_power(source)
    end
    return
end

function _get_withdrawals!(
    bus_activepower_withdrawals::Vector{Float64},
    bus_reactivepower_withdrawals::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System,
)
    # FIXME properly handle SwitchedAdmittance components
    for l in PSY.get_components(PSY.ElectricLoad, sys)
        (isa(l, PSY.FixedAdmittance) || isa(l, PSY.SwitchedAdmittance)) && continue
        !PSY.get_available(l) && continue
        bus = PSY.get_bus(l)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        # see issue #1463 in PSY
        if !isa(l, PSY.SynchronousCondenser)
            bus_activepower_withdrawals[bus_ix] += get_total_p(l)
        end
        bus_reactivepower_withdrawals[bus_ix] += get_total_q(l)
    end
    return
end

# TODO: Might need changes if we have SwitchedAdmittances
function _get_reactive_power_bound!(
    bus_reactivepower_bounds::Vector{Vector{Float64}},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System)
    for source in PSY.get_components(PSY.StaticInjection, sys)
        isa(source, PSY.ElectricLoad) && continue
        !PSY.get_available(source) && continue
        bus = PSY.get_bus(source)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        reactive_power_limits = get_reactive_power_limits_for_power_flow(source)
        if reactive_power_limits !== nothing
            bus_reactivepower_bounds[bus_ix][1] += min(0, reactive_power_limits.min)
            bus_reactivepower_bounds[bus_ix][2] += max(0, reactive_power_limits.max)
        else
            @warn("Reactive Power Bounds at Bus $(PSY.get_name(bus)) set to (-Inf, Inf)")
            bus_reactivepower_bounds[bus_ix][1] = -Inf
            bus_reactivepower_bounds[bus_ix][2] = Inf
        end
    end
end

function _initialize_bus_data!(
    bus_type::Vector{PSY.ACBusTypes},
    bus_angles::Vector{Float64},
    bus_magnitude::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System,
    correct_bustypes::Bool = false,
)
    forced_PV = must_be_PV(sys)
    possible_PV = can_be_PV(sys)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )
    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.ACBus, sys, bus_name)
        bt = PSY.get_bustype(bus)
        if (bt == PSY.ACBusTypes.PV || bt == PSY.ACBusTypes.REF) && !(bus_no in possible_PV)
            if correct_bustypes
                @info "Bus $bus_name (number $bus_no) changed from PV to PQ: no available " *
                      "sources at that bus." maxlog = PF_MAX_LOG
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
        bus_type[ix] = bt
        bus_angles[ix] = PSY.get_angle(bus)
        bus_vm = PSY.get_magnitude(bus)
        # prevent unfeasible starting values for voltage magnitude at PQ buses (for PV and REF buses we cannot do this):
        if bt == PSY.ACBusTypes.PQ && bus_vm < BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN
            @warn(
                "Initial bus voltage magnitude of $bus_vm p.u. at PQ bus $bus_name is below the plausible minimum cut-off value of $BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN p.u. and has been set to $BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN p.u.",
                maxlog = PF_MAX_LOG,
            )
            bus_vm = BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN
        elseif bt == PSY.ACBusTypes.PQ && bus_vm > BUS_VOLTAGE_MAGNITUDE_CUTOFF_MAX
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
    n_buses,
    n_branches,
    bus_lookup,
    branch_lookup,
    valid_ix,
    converged,
    loss_factors,
    correct_bustypes,
    calculate_loss_factors,
    generator_slack_participation_factors,
    voltage_stability_factors,
    calculate_voltage_stability_factors,
)
    branch_type = Vector{DataType}(undef, length(branch_lookup))
    for (ix, b) in enumerate(PNM.get_ac_branches(sys))
        branch_type[ix] = typeof(b)
    end
    bus_reactivepower_bounds = Vector{Vector{Float64}}(undef, n_buses)
    timestep_map = Dict(zip([i for i in 1:time_steps], timestep_names))
    neighbors = Vector{Set{Int}}()
    return make_powerflowdata(
        sys,
        time_steps,
        power_network_matrix,
        aux_network_matrix,
        n_buses,
        n_branches,
        bus_lookup,
        branch_lookup,
        branch_type,
        timestep_map,
        valid_ix,
        neighbors,
        converged,
        loss_factors,
        correct_bustypes,
        calculate_loss_factors,
        generator_slack_participation_factors,
        voltage_stability_factors,
        calculate_voltage_stability_factors,
    )
end

function _add_gspf_to_ijv!(
    I::Vector{Int},
    J::Vector{Int},
    V::Vector{Float64},
    sys::System,
    gsp_factors::Dict{Tuple{DataType, String}, Float64},
    bus_lookup::Dict{Int, Int},
    time_steps_iter::AbstractVector{<:Integer},
)
    for ((gen_type, gen_name), val) in gsp_factors
        val == 0.0 && continue
        gen = PSY.get_component(gen_type, sys, gen_name)
        isnothing(gen) && throw(ArgumentError("$gen_type $gen_name not found"))
        PSY.get_available(gen) || continue
        bus = PSY.get_bus(gen)
        bus_idx = bus_lookup[PSY.get_number(bus)]
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
    time_steps::Int,
    n_buses::Int,
    bus_type::Matrix{PSY.ACBusTypes},
)
    if length(generator_slack_participation_factors) == 1
        return make_bus_slack_participation_factors(
            sys,
            generator_slack_participation_factors[1],
            bus_lookup,
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
        _add_gspf_to_ijv!(I, J, V, sys, factors, bus_lookup, time_step:time_step)
    end

    bus_slack_participation_factors = sparse(I, J, V, n_buses, time_steps)
    return bus_slack_participation_factors, generator_slack_participation_factors
end

function make_bus_slack_participation_factors(
    ::System,
    ::Nothing,
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
    for gen in PSY.get_components(PSY.Generator, sys)
        if PSY.get_available(gen)
            push!(gen_buses, PSY.get_number(PSY.get_bus(gen)))
        end
    end
    # PSSe counts buses with switched shunts as PV, so we do the same here.
    for gen in PSY.get_components(PSY.SwitchedAdmittance, sys)
        if PSY.get_available(gen)
            push!(gen_buses, PSY.get_number(PSY.get_bus(gen)))
        end
    end
    return gen_buses
end

"""Return set of all bus numbers that can be PV: i.e. have an available generator,
or certain voltage regulation devices."""
function can_be_PV(sys::System)
    source_buses = must_be_PV(sys)
    for source in PSY.get_components(PSY.Source, sys)
        if PSY.get_available(source)
            push!(source_buses, PSY.get_number(PSY.get_bus(source)))
        end
    end
    return source_buses
end

function make_powerflowdata(
    sys,
    time_steps,
    power_network_matrix,
    aux_network_matrix,
    n_buses,
    n_branches,
    bus_lookup,
    branch_lookup,
    branch_type,
    timestep_map,
    valid_ix,
    neighbors,
    converged,
    loss_factors,
    calculate_loss_factors,
    correct_bustypes::Bool = false,
    generator_slack_participation_factors = nothing,
    voltage_stability_factors = nothing,
    calculate_voltage_stability_factors = nothing,
)
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = ones(Float64, n_buses)

    _initialize_bus_data!(
        bus_type,
        bus_angles,
        bus_magnitude,
        bus_lookup,
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
        sys,
    )

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    _get_withdrawals!(
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_lookup,
        sys,
    )

    # Define fields as matrices whose number of columns is equal to the number of time_steps
    bus_activepower_injection_1 = zeros(Float64, n_buses, time_steps)
    bus_reactivepower_injection_1 = zeros(Float64, n_buses, time_steps)
    bus_activepower_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    bus_reactivepower_withdrawals_1 = zeros(Float64, n_buses, time_steps)
    bus_reactivepower_bounds_1 = Matrix{Vector{Float64}}(undef, n_buses, time_steps)
    bus_magnitude_1 = ones(Float64, n_buses, time_steps)
    bus_angles_1 = zeros(Float64, n_buses, time_steps)

    # Initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= bus_magnitude
    bus_angles_1[:, 1] .= bus_angles

    bus_reactivepower_bounds = Vector{Vector{Float64}}(undef, n_buses)
    for i in 1:n_buses
        bus_reactivepower_bounds[i] = [0.0, 0.0]
    end
    _get_reactive_power_bound!(bus_reactivepower_bounds, bus_lookup, sys)
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
            time_steps,
            n_buses,
            bus_type_1,
        )

    # Initial flows are all zero
    branch_activepower_flow_from_to = zeros(Float64, n_branches, time_steps)
    branch_reactivepower_flow_from_to = zeros(Float64, n_branches, time_steps)
    branch_activepower_flow_to_from = zeros(Float64, n_branches, time_steps)
    branch_reactivepower_flow_to_from = zeros(Float64, n_branches, time_steps)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        bus_reactivepower_bounds_1,
        generator_slack_participation_factors,
        bus_slack_participation_factors,
        bus_type_1,
        branch_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_activepower_flow_from_to,
        branch_reactivepower_flow_from_to,
        branch_activepower_flow_to_from,
        branch_reactivepower_flow_to_from,
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
