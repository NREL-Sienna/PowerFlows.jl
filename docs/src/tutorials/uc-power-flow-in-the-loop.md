# Running Power Flow In The Loop with Unit Commitment

In this tutorial, you'll configure a unit commitment (UC) simulation in
[`PowerSimulations.jl`](https://nrel-sienna.github.io/PowerSimulations.jl/stable/) that calls an AC power flow solver at every dispatch interval
automatically. By the end, you'll have a simulation where the committed dispatch is
validated against the full AC network model at each hour, with results written back into
the system and optionally exported to PSS/e format.

!!! note
    This tutorial requires [`PowerSimulations.jl`](https://nrel-sienna.github.io/PowerSimulations.jl/stable/) in addition to PowerFlows.jl. Make sure
    it is installed in your Julia environment before proceeding.

If you want to run power flow independently *after* a UC simulation rather than during
it, see the tutorial [Validating a UC Dispatch with Multi-Period Power Flow](@ref).

## Setup

Load the needed packages:

```@repl uc_inloop
using PowerSystemCaseBuilder
using PowerSimulations
using PowerFlows
using PowerSystems
using HiGHS
using Logging
disable_logging(Logging.Warn)
```

Build a test system. We use the RTS-GMLC DA system, a realistic 73-bus synthetic system:

```@repl uc_inloop
sys = build_system(PSISystems, "RTS_GMLC_DA_sys")
```

!!! warning "Run the setup blocks first"
    Every code block in this tutorial shares the same REPL session (`uc_inloop`). If you
    skip or re-order blocks, variables like `sys`, `power_flow_model`, and `template_uc`
    will be undefined and you will see `UndefVarError`. Always run the setup blocks above
    before proceeding to each subsequent section.

    If any `using` statement in the setup block failed because a package was not yet
    installed, install it with `Pkg.add("PackageName")` and then **re-run the entire
    setup block** — simply installing a package does not load it into the current session.

## Configuring the Power Flow Solver

Create an [`ACPowerFlow`](@ref) solver. We attach a [`PSSEExportPowerFlow`](@ref)
exporter so that [`PowerSimulations.jl`](https://nrel-sienna.github.io/PowerSimulations.jl/stable/) automatically writes a PSS/e `.raw` file for each solved interval —
useful for downstream analysis in commercial tools. Omit the `exporter` if you don't need
PSS/e output:

```@repl uc_inloop
psse_export = PSSEExportPowerFlow(;
    psse_version = :v33,
    export_dir = joinpath("outputs", "psse"),
    overwrite = true,
)

power_flow_model = ACPowerFlow(; exporter = psse_export)
```

The `exporter` field is consumed by [`PowerSimulations.jl`](https://nrel-sienna.github.io/PowerSimulations.jl/stable/) when it calls
[`solve_and_store_power_flow!`](@ref) after each interval. The power flow solution is written
back into the `System` struct, updating generator setpoints and branch flows.

## Building the UC Problem Template

Create a `ProblemTemplate` with `power_flow_evaluation` set to the solver we just
configured. This is the key step that enables "power flow in the loop":

```@repl uc_inloop
template_uc = ProblemTemplate(
    NetworkModel(PTDFPowerModel;
        use_slacks = true,
        power_flow_evaluation = power_flow_model,
    )
)

set_device_model!(template_uc, ThermalStandard, ThermalStandardUnitCommitment)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, RenewableNonDispatch, FixedOutput)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, Line, StaticBranchUnbounded)
set_device_model!(template_uc, Transformer2W, StaticBranchUnbounded)
set_device_model!(template_uc, TapTransformer, StaticBranchUnbounded)
set_device_model!(template_uc, MonitoredLine, StaticBranch)
```

`use_slacks = true` allows the UC to remain feasible even if the AC power flow cannot
fully match the dispatch — the slack variables absorb any mismatch rather than causing
the solve to fail.

## Building and Executing the Simulation

Set up the simulation with a HiGHS optimizer and run for 3 steps (hours):

```@repl uc_inloop
optimizer = optimizer_with_attributes(HiGHS.Optimizer, "log_to_console" => false)

models = SimulationModels(;
    decision_models = [
        DecisionModel(
            template_uc,
            sys;
            name = "UC",
            optimizer = optimizer,
            system_to_file = false,
            store_variable_names = true,
        ),
    ],
)

sequence = SimulationSequence(;
    models = models,
    ini_cond_chronology = InterProblemChronology(),
)

mkpath(joinpath("outputs", "simulation"))
mkpath(joinpath("outputs", "psse"))

sim = Simulation(;
    name = "uc_with_pf",
    steps = 3,
    models = models,
    sequence = sequence,
    simulation_folder = joinpath("outputs", "simulation"),
)

build!(sim; console_level = Logging.Error, serialize = false)
```

Execute the simulation. After each UC interval solves, [`PowerSimulations.jl`](https://nrel-sienna.github.io/PowerSimulations.jl/stable/) calls
[`solve_and_store_power_flow!`](@ref) with the committed dispatch, running the AC power flow and
writing the solution back into `sys`:

!!! note
    On the first run in a fresh Julia session, this step may take several minutes while
    Julia compiles PowerSimulations.jl, JuMP, and HiGHS. Subsequent runs in the same
    session will be significantly faster.

```@repl uc_inloop
execute!(sim; enable_progress_bar = false)
```

## Inspecting Results

Load the simulation results and retrieve the UC problem results:

```@repl uc_inloop
results = SimulationResults(sim)
uc_results = get_decision_problem_results(results, "UC")
```

Read the active power output of thermal generators across the 3 solved hours:

```@repl uc_inloop
read_realized_variable(uc_results, "ActivePowerVariable__ThermalStandard")
```

Because `power_flow_evaluation` was set, the AC power flow ran at each interval and its
solution was stored in `sys`. The PSS/e `.raw` files for each solved interval are written
to the `outputs/psse/` directory, one file per interval.
