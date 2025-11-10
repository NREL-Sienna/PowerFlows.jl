function set_power_flow!(br::PSY.ACTransmission, flow::Complex)
    PSY.set_active_power_flow!(br, real(flow))
    PSY.set_reactive_power_flow!(br, imag(flow))
    return
end

function set_power_flow!(br::PNM.BranchesParallel, flow::Complex)
    for segment in br
        weight = PNM.compute_parallel_multiplier(br, PSY.get_name(segment))
        set_power_flow!(segment, flow * weight)
    end
    return
end

function set_power_flow!(br::PSY.TwoTerminalLCCLine, flow::Complex)
    PSY.set_active_power_flow!(br, real(flow))
    # TwoTerminalLCCLine does not have reactive power flow attributes (even though PFD has Q results)
    return
end

function set_power_flow!(winding::PNM.ThreeWindingTransformerWinding, flow::Complex)
    (trf, num) = (PNM.get_transformer(winding), PNM.get_winding_number(winding))
    if num == 1
        PSY.set_active_power_flow_primary!(trf, real(flow))
        PSY.set_reactive_power_flow_primary!(trf, imag(flow))
    elseif num == 2
        PSY.set_active_power_flow_secondary!(trf, real(flow))
        PSY.set_reactive_power_flow_secondary!(trf, imag(flow))
    elseif num == 3
        PSY.set_active_power_flow_tertiary!(trf, real(flow))
        PSY.set_reactive_power_flow_tertiary!(trf, imag(flow))
    else
        error("Invalid winding number: $num")
    end
    return
end

function set_voltage!(bus::PSY.ACBus, V::Complex)
    PSY.set_magnitude!(bus, abs(V))
    PSY.set_angle!(bus, angle(V))
    return
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

get_complex_voltage(bus::PSY.ACBus) = PSY.get_magnitude(bus) * exp(1im * PSY.get_angle(bus))

function get_segment_flow(
    segment::PSY.ACTransmission,
    V_from::ComplexF64,
    V_to::ComplexF64,
)
    (y11, y12, _, _) = PNM.ybus_branch_entries(segment)
    I_from = y11 * V_from + y12 * V_to
    return V_from * conj(I_from)
end
