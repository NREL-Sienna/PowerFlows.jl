precompile = @timed using PowerFlows

open("precompile_time.txt", "a") do io
    write(io, "| $(ARGS[1]) | $(precompile.time) |\n")
end

using PowerSystems
using PowerSystemCaseBuilder
using PowerFlows
using Logging
const PF = PowerFlows

configure_logging(; console_level = Logging.Info)
systems = [
    (MatpowerTestSystems, "matpower_ACTIVSg10k_sys"),
]
solvers = [PF.NewtonRaphsonACPowerFlow, PF.RobustHomotopyPowerFlow]
for (group, name) in systems
    for solver in solvers
        sys = build_system(group, name)
        try
            pf = ACPowerFlow{solver}()
            pf_data = PF.PowerFlowData(pf, sys; correct_bustypes = true)
            _, time_solve_1, _, _ = @timed PF.solve_power_flow!(pf_data; pf = pf)
            open("solve_time.txt", "a") do io
                write(
                    io,
                    "| $(ARGS[1])-$(name)-$(solver)- first solve | $(time_solve_1) |\n",
                )
            end
            pf = ACPowerFlow{solver}()
            pf_data = PF.PowerFlowData(pf, sys; correct_bustypes = true)
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
