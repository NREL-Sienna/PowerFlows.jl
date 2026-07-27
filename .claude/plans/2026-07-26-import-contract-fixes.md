# PSS/E import-contract fixes (PSCB + PowerFlowFileParser, psy6)

Upstream patches closing the five import-side gaps that break producer/consumer symmetry
with PowerFlows' object-only PSS/E exporter. Diagnosed 2026-07-26 in the PF session; see
memory note `psy6-no-component-ext-contract` for full evidence.

## Repos and branches

| Repo | Path | Branch | Commits allowed? |
|:-----|:-----|:-------|:-----------------|
| PowerSystemCaseBuilder (PSCB) | `/Users/jdlara/cache/psy6/PowerSystemCaseBuilder.jl` | `jd/import-contract-fixes` | YES — commit per task on this branch. NEVER push. |
| PowerFlowFileParser (parser) | `/Users/jdlara/cache/psy6/PowerFlowFileParser.jl` | `jd/import-contract-fixes` | YES — commit per task on this branch. NEVER push. |
| PowerFlows (PF) | `/Users/jdlara/cache/psy6/PowerFlows.jl` | `psy6-rebase` (current) | NO — leave all PF changes UNSTAGED. Never `git add`, never commit. |

PF's `test/` env has `Pkg.develop`'d the two local clones, so PF tests exercise your edits
directly (fresh `julia` process per run picks them up).

## Global constraints (bind every task)

- Commit messages: plain, descriptive, single-purpose. **NO `Co-Authored-By` trailer or any
  Claude attribution — ever.** Never push. Never touch branches other than the ones above.
- Julia style (Sienna): no `isa`/`<:` runtime checks — use dispatch; `iszero(x)` not `x == 0`;
  no ternaries; explicit `function … end` with explicit `return` for non-trivial bodies; no
  `Union{Nothing,T}` sentinel returns; terse comments (only non-obvious WHY); no dot field
  access in public-facing code — use getters.
- psy6 explicit units: PSY getters on convertible fields take a unit-system argument
  (`PSY.get_x(br, PSY.SU)`). A bare getter on a convertible field is a defect.
- Never edit files under `~/.julia/packages/` — those are read-only package checkouts. Edit
  only the three repo paths above.
- Verification per task = compile checks + the named PF ReTest subset(s). Do NOT run any
  repo's full test suite (they take an hour+); the controller does that at the end.
- Compile checks:
  - parser: `cd /Users/jdlara/cache/psy6/PowerFlowFileParser.jl && julia --project -e 'using PowerFlowFileParser'`
  - PSCB: `cd /Users/jdlara/cache/psy6/PowerSystemCaseBuilder.jl && julia --project -e 'using PowerSystemCaseBuilder'`
    (first time: `julia --project -e 'using Pkg; Pkg.instantiate()'`)
- PF acceptance runs (from `/Users/jdlara/cache/psy6/PowerFlows.jl`):
  `julia --project=test -e 'using PowerFlows; include("test/PowerFlowsTests.jl"); using .PowerFlowsTests, ReTest; retest(PowerFlowsTests, r"<REGEX>")'`
  A testset passes only when the summary row shows no Fail/Error column entries.
- If a repo has `scripts/formatter/formatter_code.jl`, run it before finishing
  (`julia --project=scripts/formatter -e 'include("scripts/formatter/formatter_code.jl")'`);
  if the formatter env fails to instantiate, note it and skip.

## Task 1 — PSCB: restore the PSS/E reimport path

**Repo:** PSCB only.

The metadata-reimport parser is disabled: `src/PowerSystemCaseBuilder.jl` line ~155 has
`#include("parsers/psse_metadata_reimport.jl")`, commented out because that file calls
`system_via_power_models`, lost in the psy6 parser migration. Every
`read_system_with_metadata` in PF's suite MethodErrors on `System(::String, ::Dict)` until
this loads. The fix set below was verified empirically by injection (VSC-reparse testset
passes with exactly these pieces).

1. Uncomment the include.
2. In `src/parsers/power_models_data.jl`, directly after `make_system`, add:
   ```julia
   """
   Construct a System from a PSS/E raw file via the PowerModels dict pipeline. Thin shim
   over [`make_system`](@ref); the metadata-reimport path threads its name-formatter kwargs
   through it.
   """
   function system_via_power_models(file_path::AbstractString; kwargs...)
       return make_system(PowerFlowFileParser.PowerModelsData(file_path); kwargs...)
   end
   ```
   (PSCB `import`s PowerFlowFileParser — check the module prefix used by neighboring code
   and match it.)
3. In `src/parsers/psse_metadata_reimport.jl`, qualify the constructor extension:
   `function System(file_path::AbstractString, md::Dict; kwargs...)` →
   `function PSY.System(...)`. PSCB has `using PowerSystems`, so the unqualified form
   creates a shadow function under Julia 1.12 (warning at include time). Match the file's
   existing `PSY.`/bare usage conventions for the body.
4. Check the rest of `psse_metadata_reimport.jl` compiles against current PSCB — it was
   verified loadable by injection, so expect no other changes; if something else is broken,
   fix minimally and report it.

**Verify:** PSCB compile check passes with no `System` shadow warning; then PF acceptance
`retest(PowerFlowsTests, r"VSC built without")` → the testset
"PSSE Exporter: a VSC built without PSS/E ext metadata re-parses (v33)" passes (1/1).

**Commit** to PSCB `jd/import-contract-fixes`.

## Task 2 — parser + PSCB: VSC DCSET → per-unit, `rated_dc_voltage`

**Repos:** parser (main edit) + PSCB (one kwarg).

PSY documents `TwoTerminalVSCLine.dc_setpoint_from/to` as per-unit — of `rated_dc_voltage`
when that side's `dc_control` is a DC-voltage mode (PSS/E converter TYPE = 1), of system
base (`baseMVA`) when it is `DC_POWER` (TYPE = 2). The parser stores raw kV/MW and never
emits `rated_dc_voltage`, so PF's DC-network lowering sees `vdc_set = 300` "pu" and the NR
diverges (3 red asserts in PF's "Parsed VSC lowers to p.u.-sane setpoints" testset).

In `src/pm_io/psse.jl`, VSC-line section (search `sub_data["dc_setpoint_from"]`, ~line 1861):

1. The block ~50 lines below already derives `base_voltage` (the voltage-controlling
   converter's DCSET, kV) and `flow_setpoint`; it computes `Zbase`, `sub_data["r"]`,
   `sub_data["pf"]`, `sub_data["if"]` from them. Reorder so `base_voltage` is derived
   BEFORE the setpoint assignment, then:
   - `sub_data["rated_dc_voltage"] = base_voltage`
   - `sub_data["dc_setpoint_from"]` = `from_bus["DCSET"] / base_voltage` when
     `from_bus["TYPE"] == 1`, else `from_bus["DCSET"] / baseMVA`; same for `to`.
   - **No sign flip.** PF's exporter inverse (`_vsc_export_dcset`) multiplies the stored pu
     value straight back by `rated_dc_voltage`/`base_power`, so magnitude-only conversion is
     the round-trip-consistent choice. Leave the existing `flow_setpoint` sign logic and the
     derived `pf`/`if` computations byte-identical.
2. The voltage-control-mode `if/elseif/else`: the `else` catches both zero AND two TYPE = 1
   converters but its message says "At least one converter … must set a voltage control".
   Split: zero voltage controllers keeps (a corrected version of) that error; two voltage
   controllers errors "exactly one converter … must control DC voltage (TYPE = 1)".
3. In PSCB `src/parsers/power_models_data.jl` `make_vscline`, add
   `rated_dc_voltage = d["rated_dc_voltage"],` to the constructor call.

**Verify:** parser + PSCB compile checks; PF acceptance
`retest(PowerFlowsTests, r"Parsed VSC lowers|round-trip preserves a VSC DC line")` — the
"Parsed VSC lowers…" testset must go fully green (all 5 asserts, including the solve), and
the VSC round-trip testset must stay green.

**Commit** each repo separately on its `jd/import-contract-fixes`.

## Task 3 — parser + PSCB: switched-shunt MODSW/SWREM onto the object

**Repos:** parser + PSCB + PF (test only, unstaged).

The parser parks `MODSW` and `SWREM`/`NREG` in `sub_data["ext"]` (psse.jl ~741-773) and
PSCB's `make_switched_shunt` (power_models_data.jl ~1839) never sets `control_mode` or
`regulated_bus_number`. Consequence: every parsed switched shunt defaults to
`control_mode = FIXED`, which PF's discrete control treats as locked — parsed shunts
silently stop regulating. PSY's `SwitchedAdmittanceControlMode` enum values ARE the MODSW
codes (PSCB already uses the analogous `TransformerControlObjective(cod)` Int-construction
pattern in `_transformer_control_fields`).

1. Parser: move MODSW out of ext into `sub_data["control_mode"] = switched_shunt["MODSW"]`.
   Add `sub_data["regulated_bus_number"]` = `SWREM` for source versions 30/32/33. For v35,
   check `pti.jl`'s v35 switched-shunt field list for the SWREG/NREG names — SWREG replaces
   SWREM as the regulated bus; NREG is the node and stays in ext. Map accordingly.
2. PSCB `make_switched_shunt`: set
   `control_mode = SwitchedAdmittanceControlMode(d["control_mode"])` and
   `regulated_bus_number = d["regulated_bus_number"]` (0 means local bus and stays 0).
   Use `get` with the current defaults only if a non-PSS/E source (matpower) reaches this
   maker without the keys — check the call sites first.
3. PF (UNSTAGED, no commit): add one round-trip testset to `test/test_psse_export.jl`
   mirroring the file's existing style (e.g. the FACTS testset at ~line 514): build a
   `System` with a `SwitchedAdmittance` with `control_mode = DISCRETE_VOLTAGE` and
   `regulated_bus_number = 7`, export v33, `read_system_with_metadata`, assert both fields
   survive. Follow the file's helper patterns (`_add_simple_bus!` etc.).

**Verify:** compile checks; PF acceptance `retest(PowerFlowsTests, r"<your new testset name>")`
plus `retest(PowerFlowsTests, r"FACTS: RMPCT|VSC built without")` for no regression.

## Task 4 — parser + PSCB: FACTS FCREG/REMOT; stop seeding `reactive_power_required`

**Repos:** parser + PSCB.

The parser's FACTS section (psse.jl ~1960-2010) never carries FCREG/REMOT and maps
`RMPCT → sub_data["reactive_power_required"]`; PSCB `make_facts`
(power_models_data.jl ~1891) passes that into the PSY field that is now documented as a
**solver output**, and never sets `regulated_bus_number`. PF's FACTS round-trip testset
fails its final assert (`get_regulated_bus_number(facts2) == 7` evaluates `0 == 7`).

1. Parser: `sub_data["regulated_bus_number"]` = `REMOT` (v30/32/33) or `FCREG` (v35 —
   pti.jl already defines it, default 0). Delete the
   `sub_data["reactive_power_required"] = facts["RMPCT"]` line; leave RMPCT in ext.
2. PSCB `make_facts`: set `regulated_bus_number = d["regulated_bus_number"]`; delete the
   `reactive_power_required = d["reactive_power_required"]` kwarg (field defaults to 0.0)
   and the now-dead `d["reactive_power_required"] < 0` validation throw.

**Verify:** compile checks; PF acceptance `retest(PowerFlowsTests, r"FACTS: RMPCT")` — the
testset must go fully green including the re-parse assert at its end.

## Task 5 — PSCB + PF: RMA/RMI radians for phase-shift control objectives

**Repos:** PSCB + PF (src + test, unstaged).

PSY documents `TransformerCircuit.control_limits` as phase-angle bounds in **radians** when
`control_objective` is an active-power (phase-shift) objective; PSS/E RMA/RMI are degrees
for those CODs. PSCB `_transformer_control_fields` (power_models_data.jl ~1351) stores them
unconverted, and PF's exporter writes the stored value back raw.

1. PSCB `_transformer_control_fields`: it already reads `cod` and normalizes band order.
   After normalization, when `TransformerControlObjective(cod)` is one of the four
   active-power objectives (`ACTIVE_POWER_FLOW`, `ACTIVE_POWER_FLOW_DISABLED`,
   `ASYMMETRIC_ACTIVE_POWER_FLOW`, `ASYMMETRIC_ACTIVE_POWER_FLOW_DISABLED` — define a local
   const tuple; PSY's `_PHASE_SHIFT_OBJECTIVES` is private and its `is_phase_shifting`
   needs a built circuit), `deg2rad` both entries of `control_limits`. When the objective is
   phase-shifting and the RMI/RMA keys are absent, default to `(min = -π, max = +π)` instead
   of the tap-band `(0.9, 1.1)` — `deg2rad(0.9)` is not a meaningful phase band.
   `controlled_quantity_limits` (VMA/VMI) are NOT angle-typed; leave untouched.
2. PF `src/psse_export.jl` (UNSTAGED, no commit): the exporter writes RMA/RMI from
   `PSY.get_control_limits(circuit)` at the 2W site (`_write_2w_transformer_record3_winding1!`)
   and the 3W site (`_collect_3w_winding_data`). Add one helper next to `_circuit_cod`:
   for a circuit whose `control_objective` is in the same four-objective set, return the
   band `rad2deg`'d; otherwise return it unchanged. Do NOT key on `PSY.is_phase_shifting`
   (true for α ≠ 0 with a voltage objective, which would wrongly convert tap-ratio bounds).
   Use it at both sites.
3. PF (UNSTAGED): add a round-trip testset to `test/test_psse_export.jl`: programmatic 2W
   transformer, `control_objective = ACTIVE_POWER_FLOW`,
   `control_limits = (min = deg2rad(-30), max = deg2rad(30))`; export v33; re-parse via
   `read_system_with_metadata`; assert the reimported circuit stores radians
   (`≈ deg2rad(±30)`) and the raw file contains the degree values.

**Verify:** compile checks (PSCB and PF: `julia --project -e 'using PowerFlows'`); PF
acceptance: your new testset green, plus
`retest(PowerFlowsTests, r"RTS regression|Exporter with case16")` for no transformer
regression.

**Commit** PSCB only; PF stays unstaged.
