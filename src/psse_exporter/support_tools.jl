# Old test helper functions, soon to be deleted

function Bus_states(sys::System)
    buses = collect(PSY.get_components(PSY.Bus, sys))
    return sort(
        DataFrame(
            "bus_name" => PSY.get_name.(buses),
            "bus_number" => PSY.get_number.(buses),
            "Vm" => PSY.get_magnitude.(buses),
            "θ" => PSY.get_angle.(buses)), [:bus_number])
end

function Bus_states(d::Dict)
    (length(d) > 1) && throw(ArgumentError, "Unimplemented for this data structure")
    return sort(select(first(values(d))["bus_results"], ["bus_name", "Vm", "θ"]))
end

function Branch_states(sys::System)
    lines = collect(PSY.get_available_components(PSY.Branch, sys))
    return sort(
        DataFrame(
            "branch_name" => PSY.get_name.(lines),
            "bus_from" => PSY.get_number.(PSY.get_from.(PSY.get_arc.(lines))),
            "bus_to" => PSY.get_number.(PSY.get_to.(PSY.get_arc.(lines))),
            "r" => PSY.get_r.(lines),
            "x" => PSY.get_x.(lines),
            #"b" => get_b.(lines),
            "P_to_from" => PSY.get_active_power_flow.(lines),
            "Q_to_from" => PSY.get_reactive_power_flow.(lines),
            "rate" => PSY.get_rating.(lines),
            #"arc" => PSY.get_arc.(lines)
        ), [:bus_from, :bus_to, :branch_name])
end

function Line_states(sys::System)
    lines = collect(PSY.get_components(PSY.Line, sys))
    return sort(
        DataFrame(
            "line_name" => PSY.get_name.(lines),
            "from_bus" => PSY.get_number.(PSY.get_from.(PSY.get_arc.(lines))),
            "to_bus" => PSY.get_number.(PSY.get_to.(PSY.get_arc.(lines))),
            "r" => PSY.get_r.(lines),
            "x" => PSY.get_x.(lines),
            "b_from" => [t.from for t in PSY.get_b.(lines)],
            "b_to" => [t.to for t in PSY.get_b.(lines)],
            "active_flow" => PSY.get_active_power_flow.(lines),
            "reactive_flow" => PSY.get_reactive_power_flow.(lines),
            "rate" => PSY.get_rating.(lines),
            #"arc" => PSY.get_arc.(lines)
        ), [:from_bus, :to_bus, "line_name"])
end

# Hacky, temporary
function getter_with_default(fn, arg, default = 0.0)
    try
        return fn(arg)
    catch
        (default isa Function) && return default(arg)
        return default
    end
end

function StaticLoad_states(sys::System)
    loads = collect(PSY.get_components(PSY.StaticLoad, sys))
    return sort(
        DataFrame(
            "load_name" => PSY.get_name.(loads),
            "load_bus" => PSY.get_number.(PSY.get_bus.(loads)),
            "constant_active_power" =>
                getter_with_default.(
                    PSY.get_constant_active_power,
                    loads,
                    PSY.get_active_power,
                ),
            "constant_reactive_power" =>
                getter_with_default.(
                    PSY.get_constant_reactive_power,
                    loads,
                    PSY.get_reactive_power,
                ),
            "impedance_active_power" =>
                getter_with_default.(PSY.get_impedance_active_power, loads),
            "impedance_reactive_power" =>
                getter_with_default.(PSY.get_impedance_reactive_power, loads),
            "current_active_power" =>
                getter_with_default.(PSY.get_current_active_power, loads),
            "current_reactive_power" =>
                getter_with_default.(PSY.get_current_reactive_power, loads),
        ),
    )
end

function FixedAdmittance_states(sys::System)
    loads = collect(PSY.get_components(PSY.FixedAdmittance, sys))
    return sort(
        DataFrame(
            "load_name" => PSY.get_name.(loads),
            "load_bus" => PSY.get_number.(PSY.get_bus.(loads)),
            "Y" => PSY.get_Y.(loads),
        ), [:load_bus])
end

function ThermalStandard_states(sys::System)
    gens = collect(PSY.get_components(PSY.ThermalStandard, sys))
    return sort(
        DataFrame(
            "Generator_name" => PSY.get_name.(gens),
            "bus_number" => PSY.get_number.(PSY.get_bus.(gens)),
            "active_power" => PSY.get_active_power.(gens),
            "reactive_power" => PSY.get_reactive_power.(gens),
            "rating" => PSY.get_rating.(gens)), [:bus_number, :Generator_name])
end

function RenewableDispatch_states(sys::System)
    gens = collect(PSY.get_components(PSY.RenewableDispatch, sys))
    return sort(
        DataFrame(
            "Generator_name" => PSY.get_name.(gens),
            "bus_number" => PSY.get_number.(PSY.get_bus.(gens)),
            "active_power" => PSY.get_active_power.(gens),
            "reactive_power" => PSY.get_reactive_power.(gens),
            "rating" => PSY.get_rating.(gens)), [:bus_number, :Generator_name])
end

function Generator_states(sys::System)
    gens = collect(PSY.get_components(PSY.Generator, sys))
    return sort(
        DataFrame(
            "Generator_name" => PSY.get_name.(gens),
            "bus_number" => PSY.get_number.(PSY.get_bus.(gens)),
            #"type" => typeof.(gens),
            "active_power" => PSY.get_active_power.(gens),
            "reactive_power" => PSY.get_reactive_power.(gens),
            "rating" => PSY.get_rating.(gens),
            "dev_base" => PSY.get_base_power.(gens)), [:bus_number, :active_power])
end

function Source_states(sys::System)
    gens = collect(PSY.get_components(PSY.Source, sys))
    return sort(
        DataFrame(
            "Generator_name" => PSY.get_name.(gens),
            "bus_number" => PSY.get_number.(PSY.get_bus.(gens)),
            #"type" => typeof.(gens),
            "active_power" => PSY.get_active_power.(gens),
            "reactive_power" => PSY.get_reactive_power.(gens),
            "rating" => 100.0,
            "dev_base" => 100.0), [:bus_number, :active_power])
end

function Transformer2W_states(sys::System)
    xfmrs = collect(PSY.get_components(PSY.Transformer2W, sys))
    return sort(
        DataFrame(
            "Transformer_name" => PSY.get_name.(xfmrs),
            "from_bus" => PSY.get_number.(PSY.get_from.(PSY.get_arc.(xfmrs))),
            "to_bus" => PSY.get_number.(PSY.get_to.(PSY.get_arc.(xfmrs))),
            #"available" => get_available.(xfmrs),
            "active_power_flow" => PSY.get_active_power_flow.(xfmrs),
            "reactive_power_flow" => PSY.get_reactive_power_flow.(xfmrs),
            #"arc" => PSY.get_arc.(xfmrs),
            "r" => PSY.get_r.(xfmrs),
            "x" => PSY.get_x.(xfmrs),
            #"primary_shunt" => get_primary_shunt.(xfmrs),
            "rate" => PSY.get_rating.(xfmrs),
        ), [:from_bus, :to_bus, :Transformer_name])
end

function TapTransformer_states(sys::System)
    xfmrs = collect(PSY.get_components(PSY.TapTransformer, sys))
    return sort(
        DataFrame(
            "Transformer_name" => PSY.get_name.(xfmrs),
            "from_bus" => PSY.get_number.(PSY.get_from.(PSY.get_arc.(xfmrs))),
            "to_bus" => PSY.get_number.(PSY.get_to.(PSY.get_arc.(xfmrs))),
            #"available" => get_available.(xfmrs),
            "active_power_flow" => PSY.get_active_power_flow.(xfmrs),
            "reactive_power_flow" => PSY.get_reactive_power_flow.(xfmrs),
            #"arc" => PSY.get_arc.(xfmrs),
            "r" => PSY.get_r.(xfmrs),
            "x" => PSY.get_x.(xfmrs),
            #"primary_shunt" => get_primary_shunt.(xfmrs),
            "rate" => PSY.get_rating.(xfmrs),
        ), [:from_bus, :to_bus, :Transformer_name])
end
