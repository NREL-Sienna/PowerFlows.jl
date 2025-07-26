#=
using Logging, DataFrames, Plots, LinearAlgebra
using PowerFlows
using PowerSystems
using PowerSystemCaseBuilder
const PSY = PowerSystems
const PF = PowerFlows
const PSB = PowerSystemCaseBuilder
=#

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
    initial_residual = (L2 = NaN, L∞ = NaN)
    final_residual = (L2 = NaN, L∞ = NaN)
    for word in ["Initial", "Final"]
        residual_norm_match =
            match(Regex("$word residual size: ([\\d.e+-]+) L2, ([\\d.e+-]+) L"), log_output)
        if residual_norm_match !== nothing
            # there should be a better way...
            if word == "Initial"
                initial_residual = (L2 = parse(Float64, residual_norm_match.captures[1]),
                    L∞ = parse(Float64, residual_norm_match.captures[2]))
            else
                final_residual = (L2 = parse(Float64, residual_norm_match.captures[1]),
                    L∞ = parse(Float64, residual_norm_match.captures[2]))
            end
        end
    end
    iteration_match = match(r"solver converged after (\d+) iterations", log_output)
    if iteration_match !== nothing
        iterations = parse(Int, iteration_match.captures[1])
    end
    converged = !any(isnan.(data.bus_activepower_withdrawals[:, 1]))
    x = zeros(2 * size(data.bus_type, 1))
    PF.update_state!(x, data, 1)
    correct_soln = isapprox(x, x_solved; norm = z -> norm(z, Inf), atol = 1e-6)
    return (
        correct = (converged && correct_soln),
        iterations = iterations,
        initial_residual_L2 = initial_residual.L2,
        initial_residual_L∞ = initial_residual.L∞,
        final_residual_L2 = final_residual.L2,
        final_residual_L∞ = final_residual.L∞,
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

function get_x_solved(sys::PSY.System)
    pf = ACPowerFlow()
    pf_data = PF.PowerFlowData(pf, sys; correct_bustypes = true)
    n = size(PF.get_bus_type(pf_data), 1)
    solve_powerflow!(pf_data; pf = pf)
    x_solved = zeros(2 * n) # x_solved is the solution we get from the power flow.
    PF.update_state!(x_solved, pf_data, 1)
    @assert !any(isnan.(pf_data.bus_activepower_withdrawals[:, 1])) "unable to calculate `x_solved`"
    return x_solved
end

function generate_df(x0s::Vector{Vector{Float64}},
    methods::Vector{DataType},
    sys::PSY.System,
)
    df = DataFrame(;
        solver = String[],
        correct = Bool[],
        iterations = Int[],
        initial_residual_L2 = Float64[],
        initial_residual_L∞ = Float64[],
        final_residual_L2 = Float64[],
        final_residual_L∞ = Float64[],
    )
    x_solved = get_x_solved(sys)

    for method in methods, x0 in x0s
        pf = ACPowerFlow{method}()
        results =
            converges_from_x0(sys, pf, x0, x_solved; maxIterations = 50, print_logs = true)
        df_row = (solver = string(method), results...)
        push!(df, df_row)
    end
    return df
end

# sample usage:
function sample_usage(trials::Int = 10, interval = (-0.05, 0.05))
    @assert trials > 0 && interval[1] < interval[2]
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys")
    n = length(get_components(PSY.ACBus, sys))
    random_in_interval(a, b) = a .+ (b - a) .* rand()
    x_solved = get_x_solved(sys)
    x0s = [
        x_solved + random_in_interval(interval...) .* rand(Float64, 2 * n) for
        _ in 1:trials
    ]
    methods = [
        NewtonRaphsonACPowerFlow,
        TrustRegionACPowerFlow,
        LevenbergMarquardtACPowerFlow,
        RobustHomotopyPowerFlow,
    ]
    df = generate_df(x0s, methods, sys)
    return df
end

using Plots;
gr();

function plot_convergence_analysis(df::DataFrame)
    # Prepare data for plotting
    plot_df = copy(df)

    # Replace iterations with -1 for incorrect solutions
    plot_df.plot_iterations = ifelse.(plot_df.correct, plot_df.iterations, -1)

    # Get unique methods for color coding
    methods = unique(plot_df.solver)
    colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]
    method_colors = Dict(zip(methods, colors[1:length(methods)]))

    # Plot 1: L2 norm of initial residual vs iterations
    p1 = scatter(;
        title = "Iterations vs Initial L2 Residual",
        xlabel = "Initial L2 Residual",
        ylabel = "Iterations to Converge",
        legend = :outerbottom,
    )

    for method in methods
        method_data = filter(row -> row.solver == method, plot_df)
        scatter!(p1,
            method_data.initial_residual_L2,
            method_data.plot_iterations;
            label = method,
            color = method_colors[method],
            alpha = 0.7,
            markersize = 4,
        )
    end

    # Add vertical line at iterations = 0 to separate failed cases

    # Plot 2: L∞ norm of initial residual vs iterations  
    p2 = scatter(;
        title = "Iterations vs. Initial L∞ Residual",
        xlabel = "Initial L∞ Residual",
        ylabel = "Iterations to Converge",
        legend = false,
    )

    for method in methods
        method_data = filter(row -> row.solver == method, plot_df)
        scatter!(p2,
            method_data.initial_residual_L∞,
            method_data.plot_iterations;
            label = method,
            color = method_colors[method],
            alpha = 0.7,
            markersize = 4,
        )
    end

    # Add vertical line at iterations = 0 to separate failed cases

    # Combine plots
    combined_plot = plot(p1, p2; layout = (2, 1), size = (800, 800))
    savefig(combined_plot, "/Users/lkiernan/Downloads/my_plot.png")
    return combined_plot
end

# Usage example:
df = sample_usage()  # Generate your data
plot_convergence_analysis(df)
# TODO: if I'm plotting iterations to converge vs initial residual, using the same x0's,
# then I end up with a bunch of overlapping points and vertical lines.
# TODO: why is RobustHomotopyPowerFlow failing on 2 points, when
# NewtonRaphson is succeeding on those same points? Nonzero reference angle maybe?
