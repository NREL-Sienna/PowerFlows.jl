# Internal

```@meta
CurrentModule = PowerFlows
DocTestSetup  = quote
    using PowerFlows
end
```

## Power Flow Data Structures

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "PowerFlowData.jl",
]
```

## LCC HVDC Parameters and Utilities

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "lcc_parameters.jl",
    "lcc_utils.jl",
]
```

## Power Flow Initialization

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "initialize_powerflow_data.jl",
]
```

## AC Power Flow - Residuals

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "ac_power_flow_residual.jl",
]
```

## AC Power Flow - Jacobian

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "ac_power_flow_jacobian.jl",
]
```

## DC Power Flow

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "solve_dc_powerflow.jl",
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
    "RobustHomotopy/homotopy_powerflow.jl",
    "RobustHomotopy/homotopy_residual.jl",
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

## PSSE Export

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "psse_export.jl",
    "psse_export_data.jl",
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

## Power Flow Types

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "powerflow_types.jl",
]
```

## Power Flow Methods

```@autodocs
Modules = [PowerFlows]
Public = false
Private = true
Pages = [
    "powerflow_method.jl",
]
```