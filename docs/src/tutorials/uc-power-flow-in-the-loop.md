# Running Power Flow In The Loop with Unit Commitment

To get started, ensure you have followed the [installation instructions](@ref "Installation"). Start Julia from the command line if you haven't already.

In this tutorial, you'll configure a unit commitment (UC) simulation in
[`PowerSimulations.jl`](https://nrel-sienna.github.io/PowerSimulations.jl/stable/) that calls an AC power flow solver at every dispatch interval
automatically. By the end, you'll have a simulation where the committed dispatch is
validated against the full AC network model at each hour, with results written back into
the system and optionally exported to PSS/e format.

!!! note
    This tutorial requires [`PowerSystemCaseBuilder.jl`](https://nrel-sienna.github.io/PowerSystemCaseBuilder.jl/stable/),
    [`PowerSimulations.jl`](https://nrel-sienna.github.io/PowerSimulations.jl/stable/), and
    [`PowerSystems.jl`](https://nrel-sienna.github.io/PowerSystems.jl/stable/) in addition to `PowerFlows.jl`. Make sure
    all are installed in your Julia environment before proceeding.

If you want to run power flow independently *after* a UC simulation rather than during
it, see the tutorial [Validating a UC Dispatch with Multi-Period Power Flow](@ref).

## Setup

!!! tip "Activate the project environment first"
    If you are following this tutorial from within the cloned PowerFlows.jl repository,
    activate the local project environment before loading packages so that Julia uses the
    local version rather than any globally-installed version:
    ```julia
    import Pkg
    Pkg.activate(".")
    Pkg.instantiate()  # first time only: downloads packages listed in Manifest.toml
    ```
    `Pkg.instantiate()` is only needed the first time you activate the project (or after
    pulling changes that update `Manifest.toml`). It is safe to skip on subsequent sessions.
    If you installed PowerFlows.jl via `Pkg.add` into your global environment (i.e. you
    are **not** working from the cloned repo), you can skip both steps. If you need to
    install the packages into your global environment first, run:
    ```julia
    import Pkg
    Pkg.add(["PowerSystemCaseBuilder", "PowerSimulations", "PowerSystems", "PowerFlows", "HiGHS"])
    ```

Load the needed packages:

```@repl uc_inloop
using PowerSystemCaseBuilder
using PowerSimulations
using PowerFlows
using PowerSystems
using HiGHS
using Logging
disable_logging(Logging.Error)
```

Build a test system. We use the 5-bus Matpower DA system, a small synthetic system:

```@repl uc_inloop
sys = build_system(PSISystems, "5_bus_matpower_DA")
```

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
the solve to fail. This is the recommended setting
where small mismatches between the simplified PTDF-based UC network model and the full AC
power flow are expected. In production workflows where you want the simulation to fail
loudly if the AC power flow cannot match the committed dispatch, set `use_slacks = false`.

## Building and Executing the Simulation

Set up the simulation with a HiGHS optimizer and run for 6 steps (hours):

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
    steps = 6,
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
execute!(sim; enable_progress_bar = true)
```

## Inspecting Results

Load the simulation results and retrieve the UC problem results:

```@repl uc_inloop
results = SimulationResults(sim)
uc_results = get_decision_problem_results(results, "UC")
```

### Commitment Decisions

Check which generators were committed (on/off) at each interval:

```@repl uc_inloop
show(read_realized_variable(uc_results, "OnVariable__ThermalStandard"), allrows = true)
```

Each row is a timestep and each column is a generator. A value of `1.0` means the unit was
committed; `0.0` means it was off. The 5-bus Matpower DA system has 6 periods of time
series data, so you will see 6 rows — one per simulated hour. In this system, all
generators are committed at every interval — the load is high enough relative to the
number of generators that the optimizer needs all of them online to meet demand. There are
no low-load periods where decommitting a unit would save money, so they all stay on. This
is expected behavior for a small test system designed for correctness checking rather than
realistic dispatch scenarios.

### Active Power Dispatch

Read the active power output of each committed generator across the 6 intervals:

```@repl uc_inloop
read_realized_variable(uc_results, "ActivePowerVariable__ThermalStandard")
```

Each value is the active power output of a generator in megawatts (MW) for that interval.
For example, a value of `532.841` means that generator was dispatched to produce roughly
533 MW at that hour — the UC optimizer chose that output level to meet the system load
while minimizing total generation cost, subject to each generator's minimum/maximum
capacity limits and ramp rate constraints. Generators that were off (`OnVariable = 0`)
will show zero output. The dispatch levels here are the UC solution from the PTDF-based
network model. At each interval, the AC power flow then re-evaluated these setpoints on
the full network — any mismatch was absorbed by the slack variables configured via
`use_slacks = true`.

### Startup Events

See which generators started up at each interval:

```@repl uc_inloop
show(read_realized_variable(uc_results, "StartVariable__ThermalStandard"), allrows = true)
```

A value of `1` indicates the unit started up at that timestep. Cross-referencing this
with the `OnVariable` results confirms the commitment trajectory: a unit that is off at
hour 1 and on at hour 2 will show a startup event at hour 2. In this system, all
generators are committed from the start and stay on throughout, so all values will be
`0.0` — no unit ever transitions from off to on, meaning no startup events occur. This is expected for such a small test system.

!!! note
    You may see values of `-0.0` alongside `0.0` in this table. These are numerically
    identical — `-0.0 == 0.0` is `true` in Julia and IEEE 754 floating-point arithmetic.
    The sign bit can be set as a side effect of certain floating-point operations (e.g.,
    multiplying a small negative number that underflows to zero). Both values mean the
    generator did not start up at that interval and can be treated the same way.

### PSS/e Export Files

Because `power_flow_evaluation` was set, the AC power flow ran at each interval and its
solution was stored in `sys`. The PSS/e `.raw` files for each solved interval are written
to the `outputs/psse/` directory, one file per interval.
