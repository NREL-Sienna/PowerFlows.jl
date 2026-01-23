# Solving a Power Flow

In this tutorial, you'll solve power flows on a 5-bus test system using three different solvers and compare their results.

```@setup basic_tutorial
using PowerSystemCaseBuilder
using PowerFlows
using PowerSystems
using Logging
```

## Building a System

Create a [`System`](@extref PowerSystems.System) from [PowerSystemCaseBuilder.jl](https://github.com/NREL-Sienna/PowerSystemCaseBuilder.jl):

```@repl basic_tutorial
sys = with_logger(SimpleLogger(stderr, Logging.Error)) do
    build_system(MatpowerTestSystems, "matpower_case5_sys")
end
```

## DC Power Flow

Create a [`DCPowerFlow`](@ref) solver:

```@repl basic_tutorial
pf_dc = DCPowerFlow()
```

Solve the power flow with [`solve_power_flow`](@ref):

```@repl basic_tutorial
dc_results = with_logger(SimpleLogger(stderr, Logging.Error)) do
    solve_power_flow(pf_dc, sys)
end
```

Look at the bus results:

```@repl basic_tutorial
dc_results["1"]["bus_results"]
```

Notice that `Vm` is 1.0 for all buses, and `Q_gen` and `Q_load` are 0. This is expected for DC power flow.

Look at the line flows:

```@repl basic_tutorial
dc_results["1"]["flow_results"]
```

Notice that `Q_from_to` and `Q_to_from` are zero.

## PTDF DC Power Flow

Create a [`PTDFDCPowerFlow`](@ref) solver:

```@repl basic_tutorial
pf_ptdf = PTDFDCPowerFlow()
```

Solve the power flow:

```@repl basic_tutorial
ptdf_results = with_logger(SimpleLogger(stderr, Logging.Error)) do
    solve_power_flow(pf_ptdf, sys)
end
```

Look at the bus results:

```@repl basic_tutorial
ptdf_results["1"]["bus_results"]
```

The results match `DCPowerFlow`. For very large systems, consider [`vPTDFDCPowerFlow`](@ref) instead.

## AC Power Flow

Create an [`ACPowerFlow`](@ref) solver:

```@repl basic_tutorial
pf_ac = ACPowerFlow()
```

Solve the power flow:

```@repl basic_tutorial
ac_results = with_logger(SimpleLogger(stderr, Logging.Error)) do
    solve_power_flow(pf_ac, sys)
end
```

Look at the bus results:

```@repl basic_tutorial
ac_results["bus_results"]
```

Notice that `Vm` now varies across buses (not all 1.0), and `Q_gen` has non-zero values.

Look at the line flows:

```@repl basic_tutorial
ac_results["flow_results"]
```

Notice that `Q_from_to` and `Q_to_from` now show reactive power flows, and `P_from_to` differs from `P_to_from` due to losses.

## When AC Power Flow Fails

Let's create a system with data issues to see what failure looks like.

```@repl basic_tutorial
bad_sys = with_logger(SimpleLogger(stderr, Logging.Error)) do
    build_system(PSITestSystems, "c_sys5_re"; add_forecasts = false)
end
```

Remove some lines and increase impedance on another to create an infeasible case:

```@repl basic_tutorial
remove_component!(Line, bad_sys, "1")
```

```@repl basic_tutorial
remove_component!(Line, bad_sys, "2")
```

```@repl basic_tutorial
br = get_component(Line, bad_sys, "6")
```

```@repl basic_tutorial
set_x!(br, 20.0)
```

```@repl basic_tutorial
set_r!(br, 2.0)
```

Now try to solve:

```@repl basic_tutorial
bad_results = solve_power_flow(ACPowerFlow(), bad_sys)
```

The solver returns `nothing` when it fails to converge. You can try a more robust solver, but if the system is truly infeasible, no solver will converge:

```@repl basic_tutorial
bad_results_robust = solve_power_flow(ACPowerFlow{TrustRegionACPowerFlow}(), bad_sys)
```

Let's decrease the resistance and impedance on that line and try again:
```@repl basic_tutorial
set_x!(br, 0.03)
set_r!(br, 0.003)
```

```@repl basic_tutorial
good_results_robust = solve_power_flow(ACPowerFlow{TrustRegionACPowerFlow}(), bad_sys)
```
