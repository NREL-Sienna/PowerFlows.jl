# Example: Using CUDA Linear Solver for Power Flow Analysis
#
# This example demonstrates how to use the GPU-accelerated cuSOLVER
# linear solver as an alternative to KLU for solving power flow problems.

using PowerFlows
using PowerSystems
using PowerSystemCaseBuilder

# Optional: Load CUDA for GPU acceleration
# If CUDA is not available, the package will fall back to KLU
try
    using CUDA
    println("CUDA is available!")
    println("CUDA devices: ", CUDA.devices())
catch e
    println("CUDA not available: ", e)
    println("Will use KLU solver instead")
end

# Load a test system
println("\nLoading power system...")
sys = build_system(PSITestSystems, "c_sys14")

# Example 1: Newton-Raphson with KLU (default)
println("\n" * "="^60)
println("Example 1: Newton-Raphson with KLU Solver")
println("="^60)

pf_klu = ACPowerFlow(
    NewtonRaphsonACPowerFlow;
    linear_solver = :klu,
    tol = 1e-9,
    maxIterations = 50
)

println("Solving with KLU...")
@time results_klu = solve_powerflow(pf_klu, sys)
println("Converged: ", results_klu["converged"])

# Example 2: Newton-Raphson with cuSOLVER (GPU)
if isdefined(Main, :CUDA) && CUDA.functional()
    println("\n" * "="^60)
    println("Example 2: Newton-Raphson with cuSOLVER (GPU)")
    println("="^60)

    pf_cuda = ACPowerFlow(
        NewtonRaphsonACPowerFlow;
        linear_solver = :cusolver,
        tol = 1e-9,
        maxIterations = 50
    )

    println("Solving with cuSOLVER...")
    @time results_cuda = solve_powerflow(pf_cuda, sys)
    println("Converged: ", results_cuda["converged"])

    # Compare results
    println("\n" * "="^60)
    println("Comparing Results")
    println("="^60)
    println("Both solvers should produce the same results (within tolerance)")
    println("Results match: ", isapprox(
        results_klu["bus_results"].Vm,
        results_cuda["bus_results"].Vm,
        rtol=1e-6
    ))
else
    println("\nSkipping cuSOLVER example (CUDA not functional)")
end

# Example 3: Trust Region with cuSOLVER
if isdefined(Main, :CUDA) && CUDA.functional()
    println("\n" * "="^60)
    println("Example 3: Trust Region with cuSOLVER")
    println("="^60)

    pf_tr_cuda = ACPowerFlow(
        TrustRegionACPowerFlow;
        linear_solver = :cusolver,
        factor = 1.0,
        eta = 0.1,
        tol = 1e-9
    )

    println("Solving with Trust Region + cuSOLVER...")
    @time results_tr = solve_powerflow(pf_tr_cuda, sys)
    println("Converged: ", results_tr["converged"])
end

# Example 4: Comparing performance on larger systems
println("\n" * "="^60)
println("Example 4: Performance Comparison (if CUDA available)")
println("="^60)

if isdefined(Main, :CUDA) && CUDA.functional()
    println("For larger systems, GPU acceleration can provide significant speedups")
    println("Try testing with larger IEEE test systems or real utility datasets")
    println("\nExample:")
    println("  sys_large = build_system(PSITestSystems, \"c_sys14\")  # Or larger system")
    println("  pf_klu = ACPowerFlow(NewtonRaphsonACPowerFlow; linear_solver=:klu)")
    println("  pf_cuda = ACPowerFlow(NewtonRaphsonACPowerFlow; linear_solver=:cusolver)")
    println("  @time solve_powerflow(pf_klu, sys_large)")
    println("  @time solve_powerflow(pf_cuda, sys_large)")
else
    println("Install CUDA.jl to compare KLU vs cuSOLVER performance:")
    println("  using Pkg")
    println("  Pkg.add(\"CUDA\")")
    println("  using CUDA")
end

println("\n" * "="^60)
println("Example completed!")
println("="^60)
