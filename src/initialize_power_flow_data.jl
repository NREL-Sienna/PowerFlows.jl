"""
Sets the fields of a PowerFlowData struct to match the given System.
"""
function initialize_power_flow_data!(
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
        pf,
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
    bus_active_power_injections = zeros(Float64, n_buses)
    bus_reactive_power_injections = zeros(Float64, n_buses)
    _get_injections!(
        pf,
        bus_active_power_injections,
        bus_reactive_power_injections,
        bus_lookup,
        reverse_bus_search_map,
        sys,
    )
    data.bus_active_power_injections[:, 1] .= bus_active_power_injections
    data.bus_reactive_power_injections[:, 1] .= bus_reactive_power_injections

    # active power withdrawals, constant current and impedance withdrawals
    bus_active_power_withdrawals = zeros(Float64, n_buses)
    bus_reactive_power_withdrawals = zeros(Float64, n_buses)
    bus_active_power_constant_current_withdrawals = zeros(Float64, n_buses)
    bus_reactive_power_constant_current_withdrawals = zeros(Float64, n_buses)
    bus_active_power_constant_impedance_withdrawals = zeros(Float64, n_buses)
    bus_reactive_power_constant_impedance_withdrawals = zeros(Float64, n_buses)
    _get_withdrawals!(
        pf,
        bus_active_power_withdrawals,
        bus_reactive_power_withdrawals,
        bus_active_power_constant_current_withdrawals,
        bus_reactive_power_constant_current_withdrawals,
        bus_active_power_constant_impedance_withdrawals,
        bus_reactive_power_constant_impedance_withdrawals,
        bus_lookup,
        reverse_bus_search_map,
        sys,
    )
    data.bus_active_power_withdrawals[:, 1] .= bus_active_power_withdrawals
    data.bus_reactive_power_withdrawals[:, 1] .= bus_reactive_power_withdrawals
    data.bus_active_power_constant_current_withdrawals[:, 1] .=
        bus_active_power_constant_current_withdrawals
    data.bus_reactive_power_constant_current_withdrawals[:, 1] .=
        bus_reactive_power_constant_current_withdrawals
    data.bus_active_power_constant_impedance_withdrawals[:, 1] .=
        bus_active_power_constant_impedance_withdrawals
    data.bus_reactive_power_constant_impedance_withdrawals[:, 1] .=
        bus_reactive_power_constant_impedance_withdrawals

    # reactive power bounds
    bus_reactive_power_bounds = Vector{Tuple{Float64, Float64}}(undef, n_buses)
    for i in 1:n_buses
        bus_reactive_power_bounds[i] = (0.0, 0.0)
    end
    _get_reactive_power_bound!(
        bus_reactive_power_bounds,
        bus_lookup,
        reverse_bus_search_map,
        sys,
    )
    data.bus_reactive_power_bounds[:, 1] .= bus_reactive_power_bounds

    # bus/generator participation factors
    # remark: everything after the 3rd argument here is contained inside data.
    make_bus_slack_participation_factors!(
        data,
        sys,
        get_slack_participation_factors(pf),
        bus_lookup,
        reverse_bus_search_map,
        length(data.time_step_map),
        n_buses,
        data.bus_type,
    )
    # LCCs: initialize parameters. For DC power flow, this also writes the fixed flows to
    # data.lcc.arc_active_power_flow_from_to and data.lcc.arc_active_power_flow_to_from.
    initialize_LCCParameters!(data, sys, bus_lookup, reverse_bus_search_map)
    # TODO VSC AC power flow model goes here.
    # LCCs and VSCs, DC only: accumulate net power into bus_hvdc_net_power.
    lcc_vsc_fixed_injections!(data, sys, bus_lookup, reverse_bus_search_map)
    # generic HVDC lines: calculate fixed flows and save to generic_hvdc_flows.
    initialize_generic_hvdc_flows!(
        data,
        sys,
        reverse_bus_search_map,
    )
    # generic HVDC lines: accumulate net power into bus_hvdc_net_power.
    hvdc_fixed_injections!(
        data,
        PSY.TwoTerminalGenericHVDCLine,
        sys,
        bus_lookup,
        reverse_bus_search_map,
    )
    # ZIP Loads, DC only: convert constant current and impedance components to constant
    # powers via assuming V = 1.0 p.u.
    handle_zip_loads!(data, pf)
    return data
end
