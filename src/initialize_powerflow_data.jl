"""
Sets the fields of a PowerFlowData struct to match the given System.
"""
function initialize_powerflow_data!(
    data::PowerFlowData,
    pf::PowerFlowEvaluationModel,
    sys::System;
    correct_bustypes = false,
)
    check_unit_setting(sys)
    nrd = get_network_reduction_data(data)
    reverse_bus_search_map = PNM.get_reverse_bus_search_map(nrd)
    bus_reduction_map = PNM.get_bus_reduction_map(nrd)
    bus_lookup = get_bus_lookup(data)
    n_buses = length(bus_lookup)

    # bus types, angles, magnitudes
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = ones(Float64, n_buses)
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
    for i in 1:size(data.bus_type, 2)
        data.bus_type[:, i] .= bus_type
    end
    data.bus_angles[:, 1] .= bus_angles
    data.bus_magnitude[:, 1] .= bus_magnitude

    # active, reactive power injections, withdrawals
    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    _get_injections!(
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_lookup,
        reverse_bus_search_map,
        sys,
    )
    data.bus_activepower_injection[:, 1] .= bus_activepower_injection
    data.bus_reactivepower_injection[:, 1] .= bus_reactivepower_injection

    # active power withdrawals, constant current and impedance withdrawals
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
    data.bus_activepower_withdrawals[:, 1] .= bus_activepower_withdrawals
    data.bus_reactivepower_withdrawals[:, 1] .= bus_reactivepower_withdrawals
    data.bus_activepower_constant_current_withdrawals[:, 1] .=
        bus_activepower_constant_current_withdrawals
    data.bus_reactivepower_constant_current_withdrawals[:, 1] .=
        bus_reactivepower_constant_current_withdrawals
    data.bus_activepower_constant_impedance_withdrawals[:, 1] .=
        bus_activepower_constant_impedance_withdrawals
    data.bus_reactivepower_constant_impedance_withdrawals[:, 1] .=
        bus_reactivepower_constant_impedance_withdrawals

    # reactive power bounds
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
    data.bus_reactivepower_bounds[:, 1] .= bus_reactivepower_bounds

    # bus/generator participation factors
    # remark: everything after the 3rd argument here is contained inside data.
    make_bus_slack_participation_factors!(
        data,
        sys,
        get_slack_participation_factors(pf),
        bus_lookup,
        reverse_bus_search_map,
        length(data.timestep_map),
        n_buses,
        data.bus_type,
    )
    # LCCs: initialize parameters. Details of LCC model depends on AC vs DC powerflow.
    initialize_LCCParameters!(data, sys, bus_lookup, reverse_bus_search_map)
    return data
end
