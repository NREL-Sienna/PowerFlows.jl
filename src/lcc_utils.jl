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
    return Vm * t * sqrt(6) / π * I_dc * cos(ϕ) * sin(α) / sin(ϕ)
end

function _update_ybus_lcc!(data::PowerFlowData, time_step::Int64)
    for (i, (fb, tb)) in enumerate(data.lcc.bus_indices)
        data.lcc.rectifier.phi[i, time_step] = _calculate_ϕ_lcc(
            data.lcc.rectifier.tap[i, time_step],
            data.lcc.rectifier.thyristor_angle[i, time_step],
            data.lcc.i_dc[i, time_step],
            data.lcc.rectifier.transformer_reactance[i],
            data.bus_magnitude[fb, time_step],
        )
        data.lcc.inverter.phi[i, time_step] = _calculate_ϕ_lcc(
            data.lcc.inverter.tap[i, time_step],
            data.lcc.inverter.thyristor_angle[i, time_step],
            -data.lcc.i_dc[i, time_step],
            data.lcc.inverter.transformer_reactance[i],
            data.bus_magnitude[tb, time_step],
        )

        rectifier_admittance = _calculate_y_lcc(
            data.lcc.rectifier.tap[i, time_step],
            data.lcc.i_dc[i, time_step],
            data.bus_magnitude[fb, time_step],
            data.lcc.rectifier.phi[i, time_step],
        )
        inverter_admittance = _calculate_y_lcc(
            data.lcc.inverter.tap[i, time_step],
            data.lcc.i_dc[i, time_step],
            data.bus_magnitude[tb, time_step],
            data.lcc.inverter.phi[i, time_step],
        )
        data.lcc.branch_admittances[i] = (rectifier_admittance, inverter_admittance)
    end
    return
end

"""
Initialize the `arcs` and `bus_indices` fields of the LCCParameters structure in the PowerFlowData.
"""
function initialize_LCC_arcs_and_buses!(
    data::PowerFlowData,
    lccs::Vector{PSY.TwoTerminalLCCLine},
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
)
    lcc_arcs = PSY.get_arc.(lccs)
    # TODO error if LCCs are involved in reductions.
    nrd = get_network_reduction_data(data)
    for (i, arc) in enumerate(lcc_arcs)
        data.lcc.arcs[i] = PNM.get_arc_tuple(arc, nrd)
        data.lcc.bus_indices[i] = (
            _get_bus_ix(
                bus_lookup,
                reverse_bus_search_map,
                PSY.get_number(PSY.get_from(arc)),
            ),
            _get_bus_ix(
                bus_lookup,
                reverse_bus_search_map,
                PSY.get_number(PSY.get_to(arc)),
            ),
        )
    end
    return
end

function initialize_LCCParameters!(
    data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData},
    sys::PSY.System,
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
)
    lccs = collect(PSY.get_components(PSY.get_available, PSY.TwoTerminalLCCLine, sys))
    isempty(lccs) && return

    initialize_LCC_arcs_and_buses!(data, lccs, bus_lookup, reverse_bus_search_map)

    # use the LCCParameter's arc_activepower_flow_{from_to/to_from} field to save the 
    # power injection at the {rectifier/inverter} respectively.

    for (i, lcc_branch) in enumerate(lccs)
        data.lcc.arc_activepower_flow_from_to[i, 1] = PSY.get_active_power_flow(lcc_branch)
        # loss curve is always in natural units.
        P_dc = with_units_base(sys, PSY.UnitSystem.NATURAL_UNITS) do
            PSY.get_active_power_flow(lcc_branch)
        end
        P_dc == 0.0 && @warn "The active_power_flow field of $(PSY.get_name(lcc_branch)) " *
              " is zero. This may lead to unexpected results in DC power flow calculations."
        lcc_loss_curve = PSY.get_loss(lcc_branch)
        loss_fcn = PSY.get_function_data(lcc_loss_curve)
        loss_constant = PSY.get_input_at_zero(lcc_loss_curve)
        loss_constant_float = isnothing(loss_constant) ? 0.0 : loss_constant
        P_loss = loss_fcn(P_dc) + loss_constant_float
        P_loss > P_dc && @warn "The loss curve of LCC $(PSY.get_name(lcc_branch)) " *
              "indicates the losses are greater than the transmitted power $P_dc. " *
              "Using 0.0 as the inverter received power."
        # convert back to per unit.
        inverter_received_power = max((P_dc - P_loss) / PSY.get_base_power(sys), 0.0)
        data.lcc.arc_activepower_flow_to_from[i, 1] = inverter_received_power
    end
    return
end

function initialize_LCCParameters!(
    data::ACPowerFlowData,
    sys::PSY.System,
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
)
    lccs = collect(PSY.get_components(PSY.get_available, PSY.TwoTerminalLCCLine, sys))
    isempty(lccs) && return

    lcc_setpoint_at_rectifier = get_lcc_setpoint_at_rectifier(data)
    @assert length(lcc_setpoint_at_rectifier) == length(lccs)
    lcc_p_set = get_lcc_p_set(data)
    lcc_i_dc = get_lcc_i_dc(data)
    lcc_dc_line_resistance = get_lcc_dc_line_resistance(data)
    lcc_rectifier_tap = get_lcc_rectifier_tap(data)
    lcc_inverter_tap = get_lcc_inverter_tap(data)
    lcc_rectifier_delay_angle = get_lcc_rectifier_thyristor_angle(data)
    lcc_inverter_extinction_angle = get_lcc_inverter_thyristor_angle(data)

    lcc_rectifier_bus = get_lcc_rectifier_bus(data)
    lcc_inverter_bus = get_lcc_inverter_bus(data)
    lcc_rectifier_transformer_reactance = get_lcc_rectifier_transformer_reactance(data)
    lcc_inverter_transformer_reactance = get_lcc_inverter_transformer_reactance(data)
    lcc_rectifier_min_alpha = get_lcc_rectifier_min_thyristor_angle(data)
    lcc_inverter_min_gamma = get_lcc_inverter_min_thyristor_angle(data)

    initialize_LCC_arcs_and_buses!(data, lccs, bus_lookup, reverse_bus_search_map)

    lcc_arcs = PSY.get_arc.(lccs)

    base_power = PSY.get_base_power(sys)
    # todo: if current set point, transform into p set point
    # lcc_p_set = I_dc_A * V_dc_V / system_base_MVA

    lcc_setpoint_at_rectifier .= (PSY.get_transfer_setpoint.(lccs) .>= 0.0)
    lcc_p_set .= abs.(PSY.get_transfer_setpoint.(lccs) ./ base_power) # only one direction is supported, no reverse flow possible
    lcc_rectifier_tap[:, 1] .= PSY.get_rectifier_tap_setting.(lccs)
    lcc_inverter_tap[:, 1] .= PSY.get_inverter_tap_setting.(lccs)
    lcc_dc_line_resistance .=
        PSY.get_r.(lccs) .+ PSY.get_rectifier_rc.(lccs) .+ PSY.get_inverter_rc.(lccs)
    lcc_i_dc .=
        (-1 .+ sqrt.(1 .+ 4 .* lcc_dc_line_resistance .* lcc_p_set)) ./
        (2 .* lcc_dc_line_resistance)
    lcc_rectifier_delay_angle[:, 1] .= PSY.get_rectifier_delay_angle.(lccs)
    lcc_inverter_extinction_angle[:, 1] .= PSY.get_inverter_extinction_angle.(lccs)
    lcc_rectifier_bus .= [
        _get_bus_ix(bus_lookup, reverse_bus_search_map, x) for
        x in PSY.get_number.(PSY.get_from.(lcc_arcs))
    ]
    lcc_inverter_bus .= [
        _get_bus_ix(bus_lookup, reverse_bus_search_map, x) for
        x in PSY.get_number.(PSY.get_to.(lcc_arcs))
    ]
    lcc_rectifier_transformer_reactance .= PSY.get_rectifier_xc.(lccs)
    lcc_inverter_transformer_reactance .= PSY.get_inverter_xc.(lccs)
    lcc_rectifier_min_alpha .=
        [x.min for x in PSY.get_rectifier_delay_angle_limits.(lccs)]
    lcc_inverter_min_gamma .=
        [x.min for x in PSY.get_inverter_extinction_angle_limits.(lccs)]
    return
end
