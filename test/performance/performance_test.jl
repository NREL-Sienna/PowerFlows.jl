precompile = @timed using PowerFlows

function is_running_on_ci()
    return get(ENV, "CI", "false") == "true" || haskey(ENV, "GITHUB_ACTIONS")
end

using Dates

pushed_to_args = false
open("precompile_time.txt", "a") do io
    if length(ARGS) == 0 && !is_running_on_ci()
        pushed_to_args = true
        push!(ARGS, "Local Test at $(Dates.now())")
    end
    write(io, "| $(ARGS[1]) | $(precompile.time) |\n")
end

using PowerSystems
using PowerSystemCaseBuilder
using PowerFlows
using Logging
import PowerFlows as PF

configure_logging(; console_level = Logging.Info)
systems = [
    (MatpowerTestSystems, "matpower_ACTIVSg10k_sys"),
]
solvers = [PF.NewtonRaphsonACPowerFlow, PF.RobustHomotopyPowerFlow]
for (group, name) in systems
    for solver in solvers
        sys = build_system(group, name)
        try
            pf = ACPowerFlow{solver}(; correct_bustypes = true)
            pf_data = PF.PowerFlowData(pf, sys)
            _, time_solve_1, _, _ = @timed PF.solve_power_flow!(pf_data; pf = pf)
            open("solve_time.txt", "a") do io
                write(
                    io,
                    "| $(ARGS[1])-$(name)-$(solver)- first solve | $(time_solve_1) |\n",
                )
            end
            pf = ACPowerFlow{solver}(; correct_bustypes = true)
            pf_data = PF.PowerFlowData(pf, sys)
            _, time_solve_2, _, _ = @timed PF.solve_power_flow!(pf_data; pf = pf)
            open("solve_time.txt", "a") do io
                write(
                    io,
                    "| $(ARGS[1])-$(name)-$(solver)- second solve | $(time_solve_2) |\n",
                )
            end
        catch e
            @error exception = (e, catch_backtrace())
            open("solve_time.txt", "a") do io
                write(io, "| $(ARGS[1])-$(name)-Solve power flow | FAILED TO TEST |\n")
            end
        end
    end
end

if !is_running_on_ci()
    for file in ["precompile_time.txt", "solve_time.txt"]
        name = replace(file, "_" => " ")[begin:(end - 4)]
        println("$name:")
        for line in eachline(open(file))
            println("\t", line)
        end
    end
    pushed_to_args && pop!(ARGS)
end