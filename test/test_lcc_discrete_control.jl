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
