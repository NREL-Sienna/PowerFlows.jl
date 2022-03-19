function psse_bus_results_compare(file_name, pf_results)
    pf_result_bus = CSV.read(file_name, DataFrame)

    v_diff = []
    angle_diff = []
    for (ix, n) in enumerate(eachrow(pf_result_bus))
        push!(v_diff, pf_results["bus_results"][ix,"Vm"] - n."Voltage (pu)")
        push!(angle_diff, pf_results["bus_results"][ix, "θ"] - (n."Angle (deg)"*(π/180)))
    end
    return v_diff, angle_diff
end

function psse_gen_results_compare(file_name, system::PSY.System)
    base_power = get_base_power(system)
    pf_result_gen = CSV.read(file_name, DataFrame)
    p_diff = []
    q_diff = []
    for row in eachrow(pf_result_gen)
        if ismissing(row."Bus  Number")
            continue
        end
        name = "generator-" * "$(row."Bus  Number")" * "-" * row."Id"
        gen = get_component(ThermalStandard, system, name)
        ap = get_active_power(gen)
        rp = get_reactive_power(gen)
        push!(p_diff, ap - row."PGen (MW)"/base_power)
        push!(q_diff, rp - row."QGen (Mvar)"/base_power)
    end
    return p_diff, q_diff
end
