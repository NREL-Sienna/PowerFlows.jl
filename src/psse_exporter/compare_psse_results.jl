using Revise
using Logging

using PowerSystems
using PowerSystemCaseBuilder
using PowerFlows

PF = PowerFlows

function PSY.get_reactive_power_limits(gen::PSY.RenewableFix) 

    gen_pf = PSY.get_power_factor(gen)

    gen_q = PSY.get_max_active_power(gen)*sqrt((1/gen_pf^2)-1)

    return (min = 0.0, max=gen_q)
end

# Compare original PSY case to version of exported .raw that has been solved in PSS/E
# Not sure if we need this for testing, may be ultimately deleted. 

include("/Users/hross2/Julia/psse_exporter/src/support_tools.jl")

# Import test cases
# sys = build_system(PSIDSystems, "WECC 240 Bus")
file_dir = "/Users/hross2/Julia/twofortybus/Marenas"
sys = with_logger(SimpleLogger(Error)) do  # Suppress system loading warnings
    System(joinpath(file_dir, "system_240[32].json"))
end
set_units_base_system!(sys, PSY.IS.UnitSystem.SYSTEM_BASE)

# Load from .raw file
file_dir2 = "/Users/hross2/Julia/psse_exporter/Raw_Export/basic/2024"
sys2 = with_logger(SimpleLogger(Error)) do  # Suppress system loading warnings
    System(joinpath(file_dir2, "basic 4_solved2.raw"))
end
set_units_base_system!(sys2, PSY.IS.UnitSystem.SYSTEM_BASE)

# DC Powerflow testing
orig_results = solve_powerflow(DCPowerFlow(), sys)
old_bus_results = Bus_states(sys)
old_branch_results = Branch_states(sys)
orig_flow_results = sort!(orig_results["1"]["flow_results"], [:bus_from, :bus_to, :line_name])
orig_bus_results = orig_results["1"]["bus_results"]

psse_bus_results = Bus_states(sys2)
psse_branch_results = Branch_states(sys2)
new_results = solve_powerflow(DCPowerFlow(), sys2)
new_flows = sort!(new_results["1"]["flow_results"], [:bus_from, :bus_to, :line_name])
new_bus_results = new_results["1"]["bus_results"]

# Getter functions compare
# @show del_Vm = old_bus_results[!, [:bus_number, :Vm]] .- psse_bus_results[!, [:bus_number, :Vm]] # e-6
# @show del_θ = old_bus_results[!, [:bus_number, :θ]] .- psse_bus_results[!, [:bus_number, :θ]] # e-7
# @show del_P_flow = old_branch_results[!, [:bus_from, :bus_to, :P_to_from]] .- psse_branch_results[!, [:bus_from, :bus_to, :P_to_from]] # not the same
# @show del_Q_flow = orig_flow_results[!, [:bus_from, :bus_to, :Q_to_from]] .- psse_branch_results[!, [:bus_from, :bus_to, :Q_to_from]] # same

# DC Powerflow results compare
# @show del_Vm = orig_bus_results[!, [:bus_number, :Vm]] .- new_bus_results[!, [:bus_number, :Vm]] #e -6
# @show del_θ = orig_bus_results[!, [:bus_number, :θ]] .- new_bus_results[!, [:bus_number, :θ]] # e-15
# @show del_P_gen = orig_bus_results[!, [:bus_number, :P_gen]] .- new_bus_results[!, [:bus_number, :P_gen]] # e-6
# @show del_P_load = orig_bus_results[!, [:bus_number, :P_load]] .- new_bus_results[!, [:bus_number, :P_load]] # same
# @show del_P_net = orig_bus_results[!, [:bus_number, :P_net]] .- new_bus_results[!, [:bus_number, :P_net]] # e-6
# @show del_P_flow = orig_flow_results[!, [:bus_from, :bus_to, :P_to_from]] .- new_flows[!, [:bus_from, :bus_to, :P_to_from]] # e-12
# @show del_Q_flow = orig_flow_results[!, [:bus_from, :bus_to, :Q_to_from]] .- new_flows[!, [:bus_from, :bus_to, :Q_to_from]]  # same

# AC Powerflow testing
orig_ac_results = solve_ac_powerflow!(sys)
orig_y_bus = PowerFlows.PowerFlowData(ACPowerFlow(), sys; check_connectivity = true).power_network_matrix.data

new_ac_results = solve_ac_powerflow!(sys2)
new_y_bus = PowerFlows.PowerFlowData(ACPowerFlow(), sys2; check_connectivity = true).power_network_matrix.data

del_y_bus = findall(orig_y_bus .!= new_y_bus)

quant_del_y_bus = DataFrame(Real = Float64[], Imag = Float64[])
for i in del_y_bus
    del_r = real(orig_y_bus[i] - new_y_bus[i])
    del_i = imag(orig_y_bus[i] - new_y_bus[i])
    push!(quant_del, [del_r, del_i])
end

@show quant_del


# gen_busses = ThermalStandard_states(sys2)
# show(gen_busses, allrows=true)

# gen_busses = sort!(append!(Generator_states(sys), Source_states(sys)), [:bus_number, :active_power])
# show(gen_busses, allrows=true)

avaialabe_gens = DataFrame("gen_name" => collect(get_name.(get_components(RenewableFix, sys))))
show(avaialabe_gens, allrows = true)

print(get_component(Source, sys, "generator-4242-ND"))

