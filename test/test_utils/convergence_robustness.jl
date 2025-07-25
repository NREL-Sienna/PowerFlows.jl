"""Runs the power flow method on the system, compares the result to the provided `x_solved`,
 and returns the tuple `(correct, iterations, x)`, where `correct` is a boolean indicating
whether the power flow converged to the correct solution, `iterations` is the number of 
iterations taken, and `x` is the state vector after the power flow has been solved. 
The `kwargs` can include any keyword arguments that `solve_powerflow!` accepts, such as 
`maxIterations` and `tol`.
"""
function converges_from_x0(sys::PSY.System,
    pf::ACPowerFlow,
    x0::Vector{Float64},
    x_solved::Vector{Float64};
    null_value::Int = -1,
    print_logs::Bool = false,
    kwargs...,
)
    data = PowerFlowData(pf, sys; correct_bustypes = true)

    # Capture logs to extract iteration count
    io = IOBuffer()
    logger = ConsoleLogger(io)
    iterations = null_value # "failed to converge" stand-in.

    with_logger(logger) do
        solve_powerflow!(data; pf = pf, x0 = x0, kwargs...)
    end

    # Extract iteration count from logs
    log_output = String(take!(io))
    print_logs && println(log_output)
    iteration_match = match(r"solver converged after (\d+) iterations", log_output)
    if iteration_match !== nothing
        iterations = parse(Int, iteration_match.captures[1])
    end
    residual_norm_match =
        match(r"Final residual size: ([\d.e+-]+) L2, ([\d.e+-]+) L", log_output)
    residual_norm = (L2 = NaN, L∞ = NaN) # default to NaN if not found
    if residual_norm_match !== nothing
        residual_norm = (L2 = parse(Float64, residual_norm_match.captures[1]),
            L∞ = parse(Float64, residual_norm_match.captures[2]))
    end

    converged = !any(isnan.(data.bus_activepower_withdrawals[:, 1]))
    x = zeros(2 * size(data.bus_type, 1))
    PF.update_state!(x, data, 1)
    correct_soln = isapprox(x, x_solved; norm = z -> norm(z, Inf), atol = 1e-6)
    return (
        correct = (converged && correct_soln),
        iterations = iterations,
        residual_norm...,
    )
end

"""Splits the state vector `x` into its components based on the bus types and return a tuple of norms."""
function split_norm(x::Vector{Float64},
    bus_types::AbstractVector{PSY.ACBusTypes},
    norm = x -> norm(x, Inf),
)
    tp = PF.partition_state(x, bus_types)
    return map(v -> norm(replace(v, NaN => 0.0)), tp)
end

# sample usage:
function sample_usage()
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys")

    # get the solution for the system
    pf = ACPowerFlow()
    pf_data = PF.PowerFlowData(pf, sys; correct_bustypes = true)
    n = size(PF.get_bus_type(pf_data), 1)
    solve_powerflow!(pf_data; pf = pf)
    x_solved = zeros(2 * n) # x_solved is the solution we get from the power flow.
    PF.update_state!(x_solved, pf_data, 1)

    @assert !any(isnan.(pf_data.bus_activepower_withdrawals[:, 1])) "unable to calculate `x_solved`"

    K = 0.1
    x0 = x_solved + randn(Float64, 2 * n) * K
    results = converges_from_x0(sys,
        ACPowerFlow{RobustHomotopyPowerFlow}(),
        x0,
        x_solved;
        null_value = -2,
        print_logs = true,
        maxIterations = 50,
    )
    println("Results: $results")
end
