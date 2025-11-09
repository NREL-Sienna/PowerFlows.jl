# CUDA Linear Solver for Power Flow

PowerFlows.jl now supports GPU-accelerated linear solvers using NVIDIA's cuSOLVER library as an alternative to the default KLU solver. This can provide significant performance improvements for large-scale power flow problems when running on systems with CUDA-capable GPUs.

## Overview

The CUDA linear solver uses `cusolverSpDcsrlsvlu` from NVIDIA's cuSOLVER library to perform sparse LU factorization and solve on the GPU. This is particularly beneficial for:

- Large power systems with thousands of buses
- Repeated power flow calculations (e.g., in time series simulations)
- Systems where GPU acceleration is available
- Both AC and DC power flow computations

## Requirements

To use the CUDA linear solver, you need:

1. A CUDA-capable NVIDIA GPU
2. CUDA Toolkit installed
3. Julia package CUDA.jl installed

## Installation

The CUDA dependency is included in PowerFlows.jl. Simply ensure CUDA.jl is available:

```julia
using Pkg
Pkg.add("CUDA")
```

## Usage

### AC Power Flow

To use the CUDA solver for AC power flow, specify `linear_solver = :cusolver` when creating an `ACPowerFlow` object:

```julia
using PowerFlows
using PowerSystems

# Load your system
sys = System("path/to/system.json")

# Create AC power flow with CUDA solver
pf = ACPowerFlow(
    NewtonRaphsonACPowerFlow;
    linear_solver = :cusolver
)

# Solve power flow
results = solve_powerflow(pf, sys)
```

### AC Power Flow with Trust Region Method

```julia
pf = ACPowerFlow(
    TrustRegionACPowerFlow;
    linear_solver = :cusolver,
    factor = 1.0,
    eta = 0.1
)

results = solve_powerflow(pf, sys)
```

### DC Power Flow

The CUDA solver also works with all DC power flow methods:

```julia
# Standard DC Power Flow
pf_dc = DCPowerFlow(linear_solver = :cusolver)
results = solve_powerflow(pf_dc, sys)

# PTDF-based DC Power Flow
pf_ptdf = PTDFDCPowerFlow(linear_solver = :cusolver)
results = solve_powerflow(pf_ptdf, sys)

# Virtual PTDF-based DC Power Flow
pf_vptdf = vPTDFDCPowerFlow(linear_solver = :cusolver)
results = solve_powerflow(pf_vptdf, sys)
```

### Comparing Solvers

You can compare the performance of KLU vs CUDA solver:

```julia
using BenchmarkTools

# KLU solver (default)
pf_klu = ACPowerFlow(NewtonRaphsonACPowerFlow; linear_solver = :klu)
@btime solve_powerflow(pf_klu, sys)

# CUDA solver
pf_cuda = ACPowerFlow(NewtonRaphsonACPowerFlow; linear_solver = :cusolver)
@btime solve_powerflow(pf_cuda, sys)
```

## Configuration Options

The `ACPowerFlow` constructor accepts a `linear_solver` parameter:

- `:klu` (default) - Uses KLU sparse LU factorization on CPU
- `:cusolver` - Uses cuSOLVER sparse LU factorization on GPU

All other power flow parameters work identically with both solvers:

```julia
ACPowerFlow(
    solver_type;
    linear_solver::Symbol = :klu,           # :klu or :cusolver
    check_reactive_power_limits::Bool = false,
    calculate_loss_factors::Bool = false,
    calculate_voltage_stability_factors::Bool = false,
    enhanced_flat_start::Bool = true,
    robust_power_flow::Bool = false,
    skip_redistribution::Bool = false
)
```

## Performance Considerations

### When to Use CUDA Solver

The CUDA solver typically provides benefits for:

- **Large systems**: Systems with >1000 buses tend to benefit more from GPU acceleration
- **Sparse matrices**: The cuSOLVER implementation is optimized for sparse matrices typical in power systems
- **Multiple solves**: Data transfer overhead is amortized across multiple iterations

### When to Use KLU Solver

The KLU solver may be preferable for:

- **Small systems**: Overhead of GPU data transfer may outweigh computation benefits
- **Systems without CUDA**: When CUDA is not available or a GPU is not present
- **Dense Jacobians**: Very dense Jacobian matrices may not benefit as much from sparse GPU solvers

## Implementation Details

### Memory Management

The CUDA solver manages GPU memory automatically:

- Matrices are transferred to GPU memory during initialization
- Factorizations are performed on the GPU
- Solutions are transferred back to CPU
- GPU resources are freed when the solver cache is garbage collected

### Numerical Stability

The CUDA solver includes the same numerical stability features as KLU:

- Iterative refinement for ill-conditioned systems
- Singular matrix detection and handling
- Pattern checking for symbolic refactorization

### Supported Methods

The CUDA linear solver is supported for:

**AC Power Flow:**
- ✅ Newton-Raphson AC Power Flow
- ✅ Trust Region AC Power Flow
- ❌ Levenberg-Marquardt (uses augmented least squares, different structure)
- ❌ Robust Homotopy (uses Hessian solvers, different approach)

**DC Power Flow:**
- ✅ Standard DC Power Flow (DCPowerFlow)
- ✅ PTDF-based DC Power Flow (PTDFDCPowerFlow)
- ✅ Virtual PTDF-based DC Power Flow (vPTDFDCPowerFlow)

## Troubleshooting

### CUDA Not Available

If you see warnings about CUDA not being available:

```julia
CUDA not available. CUSOLVERLinSolveCache will not be functional.
```

This means either:
1. CUDA.jl is not installed - run `Pkg.add("CUDA")`
2. No CUDA-capable GPU is detected
3. CUDA drivers are not properly installed

The package will still work with the default KLU solver.

### GPU Memory Issues

For very large systems, you may encounter GPU memory limitations. Monitor GPU memory usage:

```julia
using CUDA
CUDA.memory_status()
```

If memory is insufficient, either:
- Use a GPU with more memory
- Fall back to the KLU solver
- Reduce problem size if possible

## Example: Complete Workflow

```julia
using PowerFlows
using PowerSystems
using PowerSystemCaseBuilder

# Load a test system
sys = build_system(PSITestSystems, "c_sys14")

# Configure power flow with CUDA solver
pf = ACPowerFlow(
    NewtonRaphsonACPowerFlow;
    linear_solver = :cusolver,
    calculate_loss_factors = true,
    tol = 1e-9,
    maxIterations = 50
)

# Solve
results = solve_powerflow(pf, sys)

if results.converged
    println("Power flow converged!")
    println("Final bus voltages:")
    display(results.bus_results)
else
    println("Power flow did not converge")
end
```

## References

- [cuSOLVER Documentation](https://docs.nvidia.com/cuda/cusolver/)
- [CUDA.jl Documentation](https://cuda.juliagpu.org/stable/)
- [PowerFlows.jl Documentation](https://nrel-sienna.github.io/PowerFlows.jl/)
