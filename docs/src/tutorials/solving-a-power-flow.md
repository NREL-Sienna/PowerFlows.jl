# Solving a Power Flow

To get started, ensure you have followed the [installation instructions](@ref "Installation"). Start Julia from the command line if you haven't already.

In this tutorial, you'll solve power flows on a 5-bus test system using three different
solvers and compare their results.

## Building a System
To get started, load the needed packages. We're using a standard test system and want to 
keep output clean, so we adjust the logging settings to filter out a few precautionary warnings.

```@repl basic_tutorial
using PowerSystemCaseBuilder
using PowerFlows
using PowerSystems
using Logging
disable_logging(Logging.Warn)
```

Create a [`System`](@extref PowerSystems.System) from [PowerSystemCaseBuilder.jl](https://github.com/NREL-Sienna/PowerSystemCaseBuilder.jl):

```@repl basic_tutorial
sys = build_system(MatpowerTestSystems, "matpower_case5_sys")
```

!!! warning "Run the setup blocks first"
    Every code block in this tutorial shares the same REPL session (`basic_tutorial`). If
    you skip or re-order blocks, variables like `sys` and `pf_dc` will be undefined and
    you will see `UndefVarError`. Always run the setup blocks above before proceeding to
    each subsequent section.

    If any `using` statement in the setup block failed because a package was not yet
    installed, install it with `Pkg.add("PackageName")` and then **re-run the entire
    setup block** — simply installing a package does not load it into the current session.

## DC Power Flow

[`DCPowerFlow`](@ref) solves for bus voltage angles using the bus admittance matrix, 
then computes branch flows from the angle differences. Create a [`DCPowerFlow`](@ref) solver:

```@repl basic_tutorial
pf_dc = DCPowerFlow()
```

Solve the power flow with [`solve_power_flow`](@ref):

```@repl basic_tutorial
dc_results = solve_power_flow(pf_dc, sys)
```

The result is a `Dict{String, Dict{String, DataFrame}}`. The outer key is the time step
name: `"1"`. The inner dictionary stores the power flow results at that time step:
`"bus_results"` for bus data and `"flow_results"` key for AC line data.
(There's also 3rd key, `"lcc_results"`, for HVDC lines, but this sytem 
contains no such components, so the matching dataframe will be empty.) Inspect `"bus_results"`:

```@repl basic_tutorial
dc_results["1"]["bus_results"]
```

Notice that `Vm` (voltage magnitude) is 1.0 for all buses, and `Q_gen` and `Q_load` are 0.
This is expected for DC power flow, which assumes flat voltage magnitudes and ignores reactive power.

```@repl basic_tutorial
dc_results["1"]["flow_results"]
```

Likewise, `Q_from_to` and `Q_to_from` (reactive power flow on the line) are zero, for all lines.

## PTDF DC Power Flow

[`PTDFDCPowerFlow`](@ref) computes branch flows directly from bus power injections using
the Power Transfer Distribution Factor matrix, without solving for voltage angles as an
intermediate step. (This means we can omit the angle computation in contexts where we only 
care about line flows, though we don't have that option implemented here.) Create a [`PTDFDCPowerFlow`](@ref)
solver:

```@repl basic_tutorial
pf_ptdf = PTDFDCPowerFlow()
```

As before, solve the power flow with [`solve_power_flow`](@ref):

```@repl basic_tutorial
ptdf_results = solve_power_flow(pf_ptdf, sys)
```

Look at the bus results:

```@repl basic_tutorial
ptdf_results["1"]["bus_results"]
```

The results match `DCPowerFlow`, as they should: the two are mathematically equivalent. 
For very large systems where forming the full PTDF matrix would be too expensive, 
consider [`vPTDFDCPowerFlow`](@ref), which computes the same results without 
storing the dense matrix.

## AC Power Flow

Create an [`ACPowerFlow`](@ref) solver:

```@repl basic_tutorial
pf_ac = ACPowerFlow()
```

Solve the power flow:

```@repl basic_tutorial
ac_results = solve_power_flow(pf_ac, sys)
```

AC results are returned as a flat `Dict{String, DataFrame}`, with the same keys as 
before: `"bus_results"`, `"flow_results"` (AC lines), and `"lcc_results"` (HVDC lines). 
(Sienna does not support multi-period AC power flows yet.) Look at the bus results:

```@repl basic_tutorial
ac_results["bus_results"]
```

Notice that `Vm` now varies across buses (not all 1.0), and `Q_gen` has non-zero values.

Look at the line flows:

```@repl basic_tutorial
ac_results["flow_results"]
```

`Q_from_to` and `Q_to_from` now show reactive power flows, and `P_from_to` differs from
`P_to_from` due to losses.

## Comparing Results

The DC and AC active power flows will differ slightly because AC power flow accounts for
resistive losses. Compare the `P_from_to` column across all three methods:

```@repl basic_tutorial
dc_results["1"]["flow_results"][!, [:flow_name, :P_from_to]]
```

```@repl basic_tutorial
ptdf_results["1"]["flow_results"][!, [:flow_name, :P_from_to]]
```

```@repl basic_tutorial
ac_results["flow_results"][!, [:flow_name, :P_from_to]]
```

DC and PTDF-DC are identical (they are mathematically equivalent). AC values differ
because the Newton-Raphson solver finds the physically exact solution including losses.

## When AC Power Flow Fails

AC power flow is iterative and not guaranteed to converge. If it fails, `solve_power_flow`
returns `missing` and logs an error. Switch to a more robust solver and retry:

```@repl basic_tutorial
pf_tr = TrustRegionACPowerFlow()
tr_results = solve_power_flow(pf_tr, sys)
```

If the Trust Region solver also fails, [`RobustHomotopyPowerFlow`](@ref) is the most
robust option and can find solutions for systems that standard Newton methods cannot.

## Next Steps

- **How-tos**: See the how-to guides for writing results to PSS/e format, running
  multi-period power flows, and working with HVDC systems.
- **Reference**: Browse the [API reference](@ref "Public API") for all available solver
  types and their keyword arguments.
- **Explanation**: Read the explanation docs to understand when to choose each solver and
  how the algorithms work.
