function set_power_flow!(br::PSY.ACTransmission, flow::Complex)
    PSY.set_active_power_flow!(br, real(flow))
    PSY.set_reactive_power_flow!(br, imag(flow))
    return
end

function set_power_flow!(tp::Tuple{PSY.ThreeWindingTransformer, Int}, flow::Complex)
    (trf, winding) = tp
    if winding == 1
        PSY.set_active_power_flow_primary!(trf, real(flow))
        PSY.set_reactive_power_flow_primary!(trf, imag(flow))
    elseif winding == 2
        PSY.set_active_power_flow_secondary!(trf, real(flow))
        PSY.set_reactive_power_flow_secondary!(trf, imag(flow))
    elseif winding == 3
        PSY.set_active_power_flow_tertiary!(trf, real(flow))
        PSY.set_reactive_power_flow_tertiary!(trf, imag(flow))
    else
        error("Invalid winding number: $winding")
    end
    return
end

function set_voltage!(bus::PSY.ACBus, V::Complex)
    PSY.set_magnitude!(bus, abs(V))
    PSY.set_angle!(bus, angle(V))
    return
end

get_complex_voltage(bus::PSY.ACBus) = PSY.get_magnitude(bus) * exp(1im * PSY.get_angle(bus))

function calculate_segment_flow!(
    segment::Union{PSY.ACTransmission, Tuple{PSY.ThreeWindingTransformer, Int}},
    V_from::ComplexF64,
    V_to::ComplexF64,
)
    (y11, y12, _, _) = PNM.ybus_branch_entries(segment)
    I_from = y11 * V_from + y12 * V_to
    S_from = V_from * conj(I_from)
    set_power_flow!(segment, S_from)
    return
end

function calculate_segment_flow!(
    segment::Set{PSY.ACTransmission},
    V_from::ComplexF64,
    V_to::ComplexF64,
)
    for br in segment
        calculate_segment_flow!(br, V_from, V_to)
    end
    return
end
