# Tests for LCC coexistence with the discrete-control continuation. The classification
# below is the contract the checkpoint in control_continuation.jl implements: every
# per-time-step Matrix field of the LCC state must be consciously placed in exactly one
# bucket, so adding converter state without extending the checkpoint fails this test
# instead of shipping a silent rollback bug.
# i_dc classification evidence: src/lcc_utils.jl:869, written only at initialization
# in initialize_LCCParameters!, not in any solve path. Classification: CHECKPOINTED
# (solver state determined from p_set at init, read-only during solves).

const LCC_CHECKPOINTED_FIELDS = Dict(
    PF.LCCConverterParameters => [:tap, :thyristor_angle],
    PF.LCCParameters => [:i_dc],
)
const LCC_DERIVED_FIELDS = Dict(
    PF.LCCConverterParameters => [:phi],
    PF.LCCParameters => Symbol[],  # branch_admittances is a Vector (not per-ts), re-derived
)
const LCC_OUTPUT_FIELDS = Dict(
    PF.LCCConverterParameters => Symbol[],
    PF.LCCParameters => [
        :arc_active_power_flow_from_to, :arc_active_power_flow_to_from,
        :arc_reactive_power_flow_from_to, :arc_reactive_power_flow_to_from,
    ],
)
const LCC_STATIC_MATRIX_FIELDS = Dict(
    PF.LCCConverterParameters => Symbol[],
    PF.LCCParameters => [:p_set],  # setpoint schedule: input, never solver-mutated
)

@testset "LCC per-time-step state classification is exhaustive" begin
    for T in (PF.LCCParameters, PF.LCCConverterParameters)
        matrix_fields = [f for f in fieldnames(T) if fieldtype(T, f) <: Matrix]
        classified = vcat(
            LCC_CHECKPOINTED_FIELDS[T], LCC_DERIVED_FIELDS[T],
            LCC_OUTPUT_FIELDS[T], LCC_STATIC_MATRIX_FIELDS[T],
        )
        @test sort(matrix_fields) == sort(classified)
        @test allunique(classified)
    end
end

@testset "continuation checkpoint restores LCC solver state" begin
    raw = joinpath(TEST_DATA_DIR, "case5_2_lcc.raw")
    sys = make_system(PFP.PowerModelsData(raw); runchecks = false)
    pf = ACPolarPowerFlow(; check_reactive_power_limits = false)
    data = PowerFlowData(pf, sys)
    ts = 1
    @test PowerFlows.get_lcc_count(data) > 0

    # Establish a consistent derived state, then snapshot.
    PowerFlows._update_ybus_lcc!(data, ts)
    ref_rt = copy(data.lcc.rectifier.tap[:, ts])
    ref_it = copy(data.lcc.inverter.tap[:, ts])
    ref_ra = copy(data.lcc.rectifier.thyristor_angle[:, ts])
    ref_ia = copy(data.lcc.inverter.thyristor_angle[:, ts])
    ref_idc = copy(data.lcc.i_dc[:, ts])
    ref_admittances = copy(data.lcc.branch_admittances)
    snap = PowerFlows._snapshot_state(data, ts)

    # Simulate a diverged trial: garbage into every LCC solver-state column and bus state.
    data.lcc.rectifier.tap[:, ts] .= 9.9
    data.lcc.inverter.tap[:, ts] .= 9.9
    data.lcc.rectifier.thyristor_angle[:, ts] .= 1.5
    data.lcc.inverter.thyristor_angle[:, ts] .= 1.5
    data.lcc.i_dc[:, ts] .= -42.0
    data.bus_magnitude[:, ts] .= 0.5
    PowerFlows._update_ybus_lcc!(data, ts)  # caches now reflect the garbage

    PowerFlows._restore_state!(data, ts, snap)

    @test data.lcc.rectifier.tap[:, ts] == ref_rt
    @test data.lcc.inverter.tap[:, ts] == ref_it
    @test data.lcc.rectifier.thyristor_angle[:, ts] == ref_ra
    @test data.lcc.inverter.thyristor_angle[:, ts] == ref_ia
    @test data.lcc.i_dc[:, ts] == ref_idc
    # Derived caches must be re-derived at the restored state, not left stale.
    @test data.lcc.branch_admittances == ref_admittances
end
