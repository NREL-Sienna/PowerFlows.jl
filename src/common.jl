function get_total_p(l::PSY.PowerLoad)
    return PSY.get_active_power(l)
end

function get_total_q(l::PSY.PowerLoad)
    return PSY.get_reactive_power(l)
end

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

function get_total_p(l::PSY.ExponentialLoad)
    return PSY.get_active_power(l)
end

function get_total_q(l::PSY.ExponentialLoad)
    return PSY.get_reactive_power(l)
end

function _get_injections!(
    bus_activepower_injection::Vector{Float64},
    bus_reactivepower_injection::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System,
)
    sources = PSY.get_components(d -> !isa(d, PSY.ElectricLoad), PSY.StaticInjection, sys)
    for source in sources
        !PSY.get_available(source) && continue
        bus = PSY.get_bus(source)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        bus_activepower_injection[bus_ix] += PSY.get_active_power(source)
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
    loads = PSY.get_components(x -> !isa(x, PSY.FixedAdmittance), PSY.ElectricLoad, sys)
    for l in loads
        !PSY.get_available(l) && continue
        bus = PSY.get_bus(l)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        bus_activepower_withdrawals[bus_ix] += get_total_p(l)
        bus_reactivepower_withdrawals[bus_ix] += get_total_q(l)
    end
    return
end

# TODO: Might need changes if we have SwitchedAdmittances
function _get_reactive_power_bound!(
    bus_reactivepower_bounds::Vector{Vector{Float64}},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System)
    sources = PSY.get_components(d -> !isa(d, PSY.ElectricLoad), PSY.StaticInjection, sys)
    for source in sources
        !PSY.get_available(source) && continue
        bus = PSY.get_bus(source)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        reactive_power_limits = PSY.get_reactive_power_limits(source)
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
    temp_bus_map::Dict{Int, String},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System,
)
    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.ACBusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
        bus_magnitude[ix] = PSY.get_magnitude(bus)
    end
end
##############################################################################
# Matrix Methods #############################################################

function my_mul_mt(
    A::PNM.VirtualPTDF,
    x::Vector{Float64},
)
    y = zeros(length(A.axes[1]))
    for i in 1:length(A.axes[1])
        name_ = A.axes[1][i]
        y[i] = LinearAlgebra.dot(A[name_, :], x)
    end
    return y
end

function make_dc_powerflowdata(
    sys,
    time_steps,
    timestep_names,
    power_network_matrix,
    aux_network_matrix,
    extra_data,
    n_buses,
    n_branches,
    bus_lookup,
    branch_lookup,
    temp_bus_map,
    valid_ix,
)
    branch_types = Vector{DataType}(undef, length(branch_lookup))
    for (ix, b) in enumerate(PNM.get_ac_branches(sys))
        branch_types[ix] = typeof(b)
    end
    bus_reactivepower_bounds = Vector{Vector{Float64}}()
    timestep_map = Dict(zip([i for i in 1:time_steps], timestep_names))
    neighbors = Vector{Set{Int}}()
    return make_powerflowdata(
        sys,
        time_steps,
        power_network_matrix,
        aux_network_matrix,
        extra_data,
        n_buses,
        n_branches,
        bus_lookup,
        branch_lookup,
        temp_bus_map,
        branch_types,
        bus_reactivepower_bounds,
        timestep_map,
        valid_ix,
        neighbors,
    )
end

function make_powerflowdata(
    sys,
    time_steps,
    power_network_matrix,
    aux_network_matrix,
    extra_data,
    n_buses,
    n_branches,
    bus_lookup,
    branch_lookup,
    temp_bus_map,
    branch_types,
    bus_reactivepower_bounds,
    timestep_map,
    valid_ix,
    neighbors,
)
    # TODO: bus_type might need to also be a Matrix since the type can change for a particular scenario
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = zeros(Float64, n_buses)

    _initialize_bus_data!(
        bus_type,
        bus_angles,
        bus_magnitude,
        temp_bus_map,
        bus_lookup,
        sys,
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

    # initialize data
    init_1 = zeros(n_buses, time_steps)
    init_2 = zeros(n_branches, time_steps)

    # define fields as matrices whose number of columns is eqault to the number of time_steps
    bus_activepower_injection_1 = deepcopy(init_1)
    bus_reactivepower_injection_1 = deepcopy(init_1)
    bus_activepower_withdrawals_1 = deepcopy(init_1)
    bus_reactivepower_withdrawals_1 = deepcopy(init_1)
    bus_magnitude_1 = deepcopy(init_1)
    bus_angles_1 = deepcopy(init_1)
    branch_flow_values_1 = deepcopy(init_2)

    # initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= bus_magnitude
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        bus_reactivepower_bounds,
        bus_type,
        branch_types,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        timestep_map,
        valid_ix,
        power_network_matrix,
        aux_network_matrix,
        neighbors,
        extra_data,
    )
end
