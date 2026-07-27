# PSCB: bring testing to green

Grounding run (2026-07-26, PSCB clone @ 2469aba, `julia --project=test test/runtests.jl`):
**1682 pass, 0 fail, 4 error.** All four errors are the same root cause and none is a code
defect in PSCB or the parser.

## Diagnosis

The four errors:

| Testset | KeyError |
|:--|:--|
| test_parsingtestsystems Serialization/De-Serialization | `"control_mode"` |
| test_psytestsystems Serialization/De-Serialization | `"control_mode"` |
| test_transformer_parsing (case14_with_pst3w) | `"rated_dc_voltage"` |
| transformer base convention (hand-computed) | `"rated_dc_voltage"` |

The import-contract campaign made PSCB's makers **require** dict keys the campaign parser
emits (`control_mode`, `regulated_bus_number` from Task 3/4; `rated_dc_voltage` from Task 2)
— a deliberate coupled contract, blessed unguarded so missing data fails loudly. But PSCB's
`test/Project.toml` `[sources]` pins PowerFlowFileParser to the **git** psy6 rev, so the
suite ran new PSCB against the old parser. New PSCB + old parser = exactly these KeyErrors.

Secondary finding: PSCB's `test/Project.toml` sources point at `NREL-Sienna` org URLs
(IS, PowerFlowFileParser, PSY) while the package `Project.toml` points at
`Sienna-Platform` for the parser. Two different remotes for the same dep across the two
envs — works only while the branches happen to match.

## Tasks

### Task 1 — wire PSCB's test env to the local parser and re-run (the actual green gate)
- `test/Project.toml` `[sources]`: `PowerFlowFileParser = {path = "/Users/jdlara/cache/psy6/PowerFlowFileParser.jl"}`
  (mirrors what PF's test env already does). Resolve + instantiate.
- Re-run the full PSCB suite. Expected: the 4 KeyErrors disappear (the new parser emits all
  three keys); watch for any *value*-level failures newly exposed downstream of the
  serialization round-trips (none expected — 0 failures today means no assertion disagrees
  with the campaign's parse output, e.g. the RMA/RMI radians change is already live).
- Any residual failure gets diagnosed individually before any fix (no defensive `get`
  defaults in the makers — that would silently mask missing data and undo the campaign's
  loud-failure design).

### Task 2 — align the org/remote pins
- Decide the canonical org per dep (package `Project.toml` says Sienna-Platform for the
  parser; the psy6-rebase Manifest in PF also resolves Sienna-Platform). Make
  `test/Project.toml` and `Project.toml` agree for every `[sources]` entry.
- This is what upstream CI will actually exercise; the local path pin from Task 1 is
  dev-only and must be reverted in the same commit that documents the ordering (below).

### Task 3 — upstream PR ordering (why CI stays green after merge)
- The parser PR (`jd/import-contract-fixes`, adf5cb1..fa13d1a) must merge to the parser's
  psy6 branch **before** the PSCB PR (3bf4fae..2469aba). Both PSCB `[sources]` pins are
  branch pins (`rev = "psy6"`), so once the parser branch carries the commits, PSCB CI picks
  them up with no further change.
- PSCB PR description carries the contract note: makers now require `control_mode`,
  `regulated_bus_number`, `rated_dc_voltage` in the PowerModels dicts; parser ≥ the paired
  psy6 commit is a hard requirement.

### Task 4 — PSCB hygiene items deferred from the campaign (fold into the same PR or a
follow-up; none blocks green)
1. `IS.DataFormatError` sweep: ~8 pre-existing bare `DataFormatError` call sites in
   `power_system_table_data.jl` are unresolved names (UndefVarError if ever hit). Mechanical
   qualification sweep + one throw-path test.
2. Delete dead `_default_switched_shunt_v35` / `"SWREG"` table in the parser's pti.jl (noted
   in the campaign ledger; it's a foot-gun for the next reader). Parser repo, not PSCB.
3. Transformer ext duplication (final-review I2): the parser still copies COD1/CONT1/
   RMA1/RMI1/NTP1/WINDV/NOMV/MAG/R/X into `sub_data["ext"]` and PSCB attaches it verbatim —
   with RMA1/RMI1 now DEGREES in ext vs RADIANS on the circuit. Strip the keys that have
   first-class fields from the ext payload (parser side), keep genuinely unmodeled ones
   (VECGRP-class). Needs one round-trip re-check in PF afterward.
4. MODSW `UNDEFINED` round-trip flip (final-review M3): PF exports blank; pti defaults blank
   to 1 (DISCRETE_VOLTAGE) on reimport, turning a non-controlling shunt into a controlling
   one. Decide: export `0` for UNDEFINED (locked, conservative) or leave documented.
5. Ask PSY for a public `is_phase_shift_objective(::TransformerControlObjective)` predicate
   (final-review M5): three private copies of the four-objective set now exist (PSY, PSCB,
   PF); a fifth objective added in PSY silently breaks two of them.

## Execution notes

- Subagent-driven (Sonnet implementers, per-task review), same rules as the campaign:
  commits only on `jd/import-contract-fixes` in the two clones, PF stays unstaged, no
  pushes, no attribution trailers.
- Full PSCB suite is the gate for Tasks 1-2 (it's ~15 min, acceptable as the acceptance run
  at controller level, background). PF's export testsets re-run after Task 4.3 only.
