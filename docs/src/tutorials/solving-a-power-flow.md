# Solving a Power Flow

```@setup basic_tutorial
using PowerSystemCaseBuilder
using PowerFlows
using PowerSystems
```

## Step 1: Build a System

Create a [System](@extref System) object. Here, we'll use a pre-made case from `PowerSystemCaseBuilder.jl`.

```@repl basic_tutorial
sys = build_system(MatpowerTestSystems, "matpower_case5_sys")
```

`display(sys)` reveals that this system has 5 [thermal generators](@extref ThermalStandard) and 3 [fixed loads](@extref PowerLoad). We could go check their names, locations, set points with the various getters...but for that, refer to the PowerSystems.jl tutorial. 

## Step 2: Define the Power Flow
Initialize a power flow solver. Here, we'll highlight two options, [ACPowerFlow](@ref) and [DCPowerFlow](@ref)
```@repl basic_tutorial
pf_ac = ACPowerFlow()
pf_dc = DCPowerFlow()
```
There's also [PTDFDCPowerFlow](@ref) and [vPTDFDCPowerFlow](@ref).

## Step 3: Solve

```@repl basic_tutorial
ac_results = solve_powerflow(pf_ac, sys)
dc_results = solve_powerflow(pf_dc, sys)
```
Now, we can inspect the power injections/withdrawals at the buses:
```@repl basic_tutorial
dc_results["1"]["bus_results"]
ac_results["bus_results"]
```

Notice that the `P_load` column is exactly the same between the two results: our power flow 
solves for generator setpoints, and leaves the load numbers alone. The `P_gen` column is 
*almost* the same between the two results: the differences reflect the choice of model.

Similarly, we can compare the line flows:
```@repl basic_tutorial
dc_results["1"]["flow_results"]
ac_results["flow_results"]
```

## Step 4 (Optional): Adjust the System and Solve Again

Let's increase the load at bus 2 and then solve a power flow again. First, inspect the current load:
```@repl basic_tutorial
load_2 = get_component(PowerLoad, sys, "bus2")
```

The load's active power and its max active power are both `3` (per unit). 
Let's increase both of those numbers by 10%:
```@repl basic_tutorial
set_max_active_power!(load_2, get_max_active_power(load_2) * 1.10)
set_active_power!(load_2, get_active_power(load_2) * 1.10)
```

Now we can solve for this higher load scenario, using the same power flow solvers as before:
```@repl basic_tutorial
new_ac_results = solve_powerflow(pf_ac, sys)
new_dc_results = solve_powerflow(pf_dc, sys)
```
