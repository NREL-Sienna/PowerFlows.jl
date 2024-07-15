# NOTE this is a standalone script -- neither a part of the library codebase nor a part of the standard test suite.
# Invoke with `include("path/to/test_powerflow.jl")` in an appropriate environment.
# TODO remove when superseded by encapsulated utility functions and/or standard test suite.

using Revise
using Logging

using PowerSystems
using PowerSystemCaseBuilder
using PowerFlows

const PF = PowerFlows

# Load test system
#show_systems()
# sys = build_system(PSSEParsingTestSystems, "pti_case3_sys")

# sys = build_system(PSIDSystems, "psid_11bus_andes")
# andes_loads = collect(get_name.(get_components(StandardLoad, sys)))
# andes_load = get_component(StandardLoad, sys, "load91")

# sys = System(joinpath(file_dir, "PSCAD_VALIDATION_RAW.raw"); bus_name_formatter = x -> strip(string(x["name"])) * "-" * string(x["index"]))
# sys = build_system(PSIDSystems, "WECC 240 Bus")

file_dir = joinpath(PF.DATA_DIR, "twofortybus", "Marenas")
sys = with_logger(SimpleLogger(Error)) do  # Suppress system loading warnings
    System(joinpath(file_dir, "system_240[32].json"))
end

set_units_base_system!(sys, UnitSystem.SYSTEM_BASE)

tfy_therms = sort!(collect(get_name.(get_components(ThermalStandard, sys))))
tfy_re_gen = get_component(RenewableDispatch, sys, "generator-2911-S-gfl")
#set_units_base_system!(sys, UnitSystem.NATURAL_UNITS)
tfy_therm_gen = get_component(ThermalStandard, sys, "generator-1436-C")
#set_units_base_system!(sys, UnitSystem.DEVICE_BASE)
tfy_therm_gen = get_component(ThermalStandard, sys, "generator-1431-N")
tfy_oad = get_component(StandardLoad, sys, "load12021")

# Solve powerlfow and get results
orig_results = solve_powerflow(DCPowerFlow(), sys)
orig_flows = sort!(orig_results["1"]["flow_results"], [:bus_from, :bus_to, :line_name])
orig_buss_results = orig_results["1"]["bus_results"]

busses_before = PF.Bus_states(sys)
lines_before = PF.Line_states(sys)
loads_before = PF.StandardLoad_states(sys)
fixed_admit_before = PF.FixedAdmittance_states(sys)
therm_gens_before = PF.ThermalStandard_states(sys)
gens_before =
    sort!(
        append!(PF.Generator_states(sys), PF.Source_states(sys)),
        [:bus_number, :active_power],
    )
xfmrs_before = sort!(
    append!(PF.Transformer2W_states(sys), PF.TapTransformer_states(sys)),
    [:from_bus, :to_bus],
)
shunt_before = PF.FixedAdmittance_states(sys)

# Write to .raw file
PF.Write_Sienna2PSSE(
    sys,
    "basic",
    2024;
    export_location = joinpath(PF.DATA_DIR, "export"),
    v33 = true,
)

# Load from .raw file
file_dir2 = joinpath(PF.DATA_DIR, "export", "Raw_Export", "basic", "2024")
# sys2 = System(joinpath(file_dir2, "basic.raw"))
sys2 = with_logger(SimpleLogger(Error)) do  # Suppress system loading warnings
    System(joinpath(file_dir2, "basic.raw"))
end
set_units_base_system!(sys, UnitSystem.SYSTEM_BASE)
# Solve popwerflpw and get new results
new_results = solve_powerflow(DCPowerFlow(), sys2)
new_flows = sort!(new_results["1"]["flow_results"], [:bus_from, :bus_to, :line_name])
new_bus_results = new_results["1"]["bus_results"]

busses_after = PF.Bus_states(sys2)
lines_after = PF.Line_states(sys2)
loads_after = PF.StandardLoad_states(sys2)
therm_gens_after = PF.ThermalStandard_states(sys2)
gens_after = PF.Generator_states(sys2)
xfmrs_after = sort!(
    append!(PF.Transformer2W_states(sys2), PF.TapTransformer_states(sys2)),
    [:from_bus, :to_bus],
)
shunt_after = PF.FixedAdmittance_states(sys2)

# Compare results
# println("Orig:", orig_flows)
# println("New:", new_flows)

# println("Orig Gen:", gens_before)
# println("New Gen:", gens_after)

# println("Orig Gen:", therm_gens_before)
# println("New Gen:", therm_gens_after)

# println("Orig Xfmr:", xfmrs_before)
# println("New Xfmr:", xfmrs_after)

println("Orig Line:", lines_before)
println("New Lines:", lines_after)

# println("Orig Loads:", loads_before)
# println("New Loads:", loads_after)

# println("Orig Bus:", orig_buss_results)
# println("New Bus:", new_bus_results)

orig_compare = PF.compare_begin_to_final(orig_flows, new_flows)
bus_compare = PF.compare_begin_to_final(orig_buss_results, new_bus_results)
line_compare = PF.compare_begin_to_final(lines_before, lines_after)
load_compare = PF.compare_begin_to_final(loads_before, loads_after)
gens_compare = PF.compare_begin_to_final(gens_before, gens_after)
# therm_gens_compare = PF.compare_begin_to_final(therm_gens_before, therm_gens_after)
# xfmr_compare = PF.compare_begin_to_final(xfmrs_before, xfmrs_after)
shunt_compare = PF.compare_begin_to_final(shunt_before, shunt_after)

println("Flow Results:", orig_compare[!, 17:end]) #e-12 error
# println("Bus Results:", bus_compare[!, 17:end]) # e-15 error
# println("Lines:", line_compare[!, 12:end])  #!!! NO REAL AND REACTIVE POWER FLOW
# println("Loads:", load_compare[!, 15:end]) # 100% Match
# println("Gens:", gens_compare[!, 9:end]) #e-15 err
# println("Xfmrs:", xfmr_compare[!, 11:end]) # !!! NO REAL AND REACTIVE POWER FLOW
# println("Fixed Admit:", shunt_compare) # 100% Match

#set_units_base_system!(sys, PSY.IS.UnitSystem.SYSTEM_BASE)
# set_units_base_system!(sys, PSY.IS.UnitSystem.NATURAL_UNITS)
# @show sum(collect(get_active_power.(get_available_components(Generator, sys))))
# @show sum(collect(get_active_power.(get_available_components(Generator, sys2))))
# @show sum(collect(get_reactive_power.(get_available_components(Generator, sys))))
# @show sum(collect(get_reactive_power.(get_available_components(Generator, sys2))))
# @show sum(collect(get_constant_active_power.(get_available_components(StandardLoad, sys))))
# @show sum(collect(get_constant_active_power.(get_available_components(StandardLoad, sys2))))
# @show sum(collect(get_constant_reactive_power.(get_available_components(StandardLoad, sys))))
# @show sum(collect(get_constant_reactive_power.(get_available_components(StandardLoad, sys2))))

# @show sum(collect(get_max_active_power.(get_available_components(Generator, sys2))))
