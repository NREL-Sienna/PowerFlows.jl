# Internal

```@meta
CurrentModule = PowerFlows
DocTestSetup  = quote
    using PowerFlows
end
```

## Power Flow Types

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "power_flow_types.jl",
]
```

# PowerFlowData

## Struct and Type Definitions 

```@autodocs
Modules = [PowerFlows]
Public = true
Private = true
Pages = [
    "PowerFlowData.jl",
]
```

## Solving a PowerFlowData instance

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "solve_ac_power_flow.jl",
    "solve_dc_power_flow.jl"
]
```

## Manipulating a PowerFlowData instance

```@autodocs
Modules = [PowerFlows]
Public = true
Private = true
Pages = [
    "state_indexing_helpers.jl",
    "initialize_power_flow_data.jl",
    "power_flow_setup.jl",
]
```

# LCC HVDC Parameters and Utilities

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "lcc_parameters.jl",
    "lcc_utils.jl",
]
```


# AC Power Flow

## Residuals

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "ac_power_flow_residual.jl",
]
```

## Jacobian

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "ac_power_flow_jacobian.jl",
]
```

## Iterative Methods

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "power_flow_method.jl",
]
```

## Robust Homotopy Method

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "RobustHomotopy/homotopy_hessian.jl",
    "RobustHomotopy/homotopy_initialization.jl",
    "RobustHomotopy/homotopy_power_flow.jl",
    "RobustHomotopy/homotopy_residual.jl",
]
```

## Levenberg-Marquardt Method

```@autodocs
Modules = [PowerFlows]
Public = true
Private = true
Pages = [
    "levenberg-marquardt.jl",
]
```

# Linear Algebra Backends

## Robust Homotopy
```@autodocs
Modules = [PowerFlows]
Public = true
Private = true
Pages = [
    "RobustHomotopy/HessianSolver/cholesky_solver.jl",
    "RobustHomotopy/HessianSolver/fixed_structure_CHOLMOD.jl",
    "RobustHomotopy/HessianSolver/hessian_solver.jl",
    "RobustHomotopy/HessianSolver/KLU_hessian_solver.jl",
]
```

## Newton-Raphson
```@autodocs
Modules = [PowerFlows]
Public = true
Private = true
Pages = [
    "LinearSolverCache/klu_linear_solver.jl",
    "LinearSolverCache/linear_solver_cache.jl",
]
```

# Misc.

## PSSE Export

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "psse_export.jl",
]
```

## Post-Processing

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "post_processing.jl",
]
```

## Power Systems Utilities

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "powersystems_utils.jl",
]
```

## Common Utilities and Definitions

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "common.jl",
    "definitions.jl",
    "psi_utils.jl",
]
```

