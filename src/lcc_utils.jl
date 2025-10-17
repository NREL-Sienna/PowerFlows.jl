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
    for (i, (fb, tb)) in enumerate(zip(data.lcc.rectifier.bus, data.lcc.inverter.bus))
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
        data.lcc.branch_admittances[(fb, tb)] = (rectifier_admittance, inverter_admittance)
    end
    return
end

function initialize_LCCParameters!(
    data::PowerFlowData,
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

    lcc_arcs = PSY.get_arc.(lccs)
    for arc in lcc_arcs
        push!(
            data.lcc.arcs,
            (PSY.get_number(PSY.get_from(arc)), PSY.get_number(PSY.get_to(arc))),
        )
    end

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
