using DataFrames
PSY = PowerSystems

function Bus_states(sys::System)
    buses = collect(get_components(Bus, sys))
    return sort(
        DataFrame(
            "bus_name" => get_name.(buses),
            "bus_number" => get_number.(buses),
            "Vm" => get_magnitude.(buses),
            "θ" => get_angle.(buses)), [:bus_number])
end

function Bus_states(d::Dict)
    (length(d) > 1) && throw(ArgumentError, "Unimplemented for this data structure")
    return sort(select(first(values(d))["bus_results"], ["bus_name", "Vm", "θ"]))
end

function Branch_states(sys::System)
    lines = collect(get_available_components(Branch, sys))
    return sort(
        DataFrame(
            "branch_name" => get_name.(lines),
            "bus_from" => get_number.(get_from.(get_arc.(lines))),
            "bus_to" => get_number.(get_to.(get_arc.(lines))),
            "r" => get_r.(lines),
            "x" => get_x.(lines),
            #"b" => get_b.(lines),
            "P_to_from" => get_active_power_flow.(lines),
            "Q_to_from" => get_reactive_power_flow.(lines),
            "rate" => get_rate.(lines),
            #"arc" => get_arc.(lines)
        ), [:bus_from, :bus_to, :branch_name])
end

function Line_states(sys::System)
    lines = collect(get_components(Line, sys))
    return sort(
        DataFrame(
            "line_name" => get_name.(lines),
            "from_bus" => get_number.(get_from.(get_arc.(lines))),
            "to_bus" => get_number.(get_to.(get_arc.(lines))),
            "r" => get_r.(lines),
            "x" => get_x.(lines),
            #"b" => get_b.(lines),
            "active_flow" => get_active_power_flow.(lines),
            "reactive_flow" => get_reactive_power_flow.(lines),
            "rate" => get_rate.(lines),
            #"arc" => get_arc.(lines)
        ), [:from_bus, :to_bus, "line_name"])
end

function StandardLoad_states(sys::System)
    loads = collect(get_components(StandardLoad, sys))
    return sort(
        DataFrame(
            "load_name" => get_name.(loads),
            "load_bus" => get_number.(get_bus.(loads)),
            "constant_active_power" => get_constant_active_power.(loads),
            "constant_reactive_power" => get_constant_reactive_power.(loads),
            "impedance_active_power" => get_impedance_active_power.(loads),
            "impedance reactive_power" => get_impedance_reactive_power.(loads),
            "current_active_power" => get_current_active_power.(loads),
            "current_reactive_power" => get_current_reactive_power.(loads),
        ),
    )
end

function FixedAdmittance_states(sys::System)
    loads = collect(get_components(FixedAdmittance, sys))
    return sort(
        DataFrame(
            "load_name" => get_name.(loads),
            "load_bus" => get_number.(get_bus.(loads)),
            "Y" => get_Y.(loads),
        ), [:load_bus])
end

function ThermalStandard_states(sys::System)
    gens = collect(get_components(ThermalStandard, sys))
    return sort(
        DataFrame(
            "Generator_name" => get_name.(gens),
            "bus_number" => get_number.(get_bus.(gens)),
            "active_power" => get_active_power.(gens),
            "reactive_power" => get_reactive_power.(gens),
            "rating" => get_rating.(gens)), [:bus_number, :Generator_name])
end

function RenewableDispatch_states(sys::System)
    gens = collect(get_components(RenewableDispatch, sys))
    return sort(
        DataFrame(
            "Generator_name" => get_name.(gens),
            "bus_number" => get_number.(get_bus.(gens)),
            "active_power" => get_active_power.(gens),
            "reactive_power" => get_reactive_power.(gens),
            "rating" => get_rating.(gens)), [:bus_number, :Generator_name])
end

function Generator_states(sys::System)
    gens = collect(get_components(Generator, sys))
    return sort(
        DataFrame(
            "Generator_name" => get_name.(gens),
            "bus_number" => get_number.(get_bus.(gens)),
            #"type" => typeof.(gens),
            "active_power" => get_active_power.(gens),
            "reactive_power" => get_reactive_power.(gens),
            "rating" => get_rating.(gens),
            "dev_base" => get_base_power.(gens)), [:bus_number, :active_power])
end

function Source_states(sys::System)
    gens = collect(get_components(Source, sys))
    return sort(
        DataFrame(
            "Generator_name" => get_name.(gens),
            "bus_number" => get_number.(get_bus.(gens)),
            #"type" => typeof.(gens),
            "active_power" => get_active_power.(gens),
            "reactive_power" => get_reactive_power.(gens),
            "rating" => 100.0,
            "dev_base" => 100.0), [:bus_number, :active_power])
end

function Transformer2W_states(sys::System)
    xfmrs = collect(get_components(Transformer2W, sys))
    return sort(
        DataFrame(
            "Transformer_name" => get_name.(xfmrs),
            "from_bus" => get_number.(get_from.(get_arc.(xfmrs))),
            "to_bus" => get_number.(get_to.(get_arc.(xfmrs))),
            #"available" => get_available.(xfmrs),
            "active_power_flow" => get_active_power_flow.(xfmrs),
            "reactive_power_flow" => get_reactive_power_flow.(xfmrs),
            #"arc" => get_arc.(xfmrs),
            "r" => get_r.(xfmrs),
            "x" => get_x.(xfmrs),
            #"primary_shunt" => get_primary_shunt.(xfmrs),
            "rate" => get_rate.(xfmrs),
        ), [:from_bus, :to_bus, :Transformer_name])
end

function TapTransformer_states(sys::System)
    xfmrs = collect(get_components(TapTransformer, sys))
    return sort(
        DataFrame(
            "Transformer_name" => get_name.(xfmrs),
            "from_bus" => get_number.(get_from.(get_arc.(xfmrs))),
            "to_bus" => get_number.(get_to.(get_arc.(xfmrs))),
            #"available" => get_available.(xfmrs),
            "active_power_flow" => get_active_power_flow.(xfmrs),
            "reactive_power_flow" => get_reactive_power_flow.(xfmrs),
            #"arc" => get_arc.(xfmrs),
            "r" => get_r.(xfmrs),
            "x" => get_x.(xfmrs),
            #"primary_shunt" => get_primary_shunt.(xfmrs),
            "rate" => get_rate.(xfmrs),
        ), [:from_bus, :to_bus, :Transformer_name])
end

function compare_begin_to_final(frame_before::DataFrame, frame_after::DataFrame)
    #@assert frame_before[!, 1] == frame_after[!, 1]
    orig_names = names(frame_before)
    del_names =
        DataFrame(["Δ" * name => zeros(nrow(frame_before)) for name in orig_names[2:end]])
    frame_combined = hcat(
        rename(frame_before, [n => n * "_before" for n in names(frame_before)[2:end]]),
        select(
            rename(frame_after, [n => n * "_after" for n in names(frame_after)]),
            Not(1),
        ),
    )
    for name in orig_names[2:end]
        del_names[!, Symbol("Δ" * name)] =
            frame_combined[!, Symbol(name * "_before")] .-
            frame_combined[!, Symbol(name * "_after")]
    end
    return final_frame = hcat(frame_combined, del_names)
end

function get_PSSE_status(bool::Bool)
    if (bool)
        return 1
    else
        return 0
    end
end

function get_psse_gen_wmod(gen::GEN) where {GEN <: Union{PSY.Generator, PSY.Storage}}
    if (occursin("Wind", PSY.get_name(gen)))
        return 1
    else
        return 0
    end
end
