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

@testset "LCC + control_discrete_devices constructs at any time_steps" begin
    raw = joinpath(TEST_DATA_DIR, "case5_2_lcc.raw")
    sys = make_system(PFP.PowerModelsData(raw); runchecks = false)
    for nts in (1, 3)
        pf = ACPolarPowerFlow(; control_discrete_devices = true, time_steps = nts)
        data = PowerFlowData(pf, sys)   # must not throw
        @test PowerFlows.get_lcc_count(data) > 0
    end
end

"""Parse the bundled two-LCC fixture and add enrollable controlled devices: a stepping
switched shunt and a shunt FACTS device at bus 101 (PQ, 230 kV, largest load). The
fixture carries no transformers, so tap devices are exercised by the existing
discrete-control suites; the tap/LCC interaction is limited to the shared checkpoint,
which is device-independent. `p_set_mw` overrides both LCC transfer setpoints (0.0
exercises the i_dc = 0 tap-pinning branch).

Every branch in `case5_2_lcc.raw` carries x = 1e-4 pu (a minimal synthetic topology sized
for LCC solving, not for realistic bus-to-bus electrical distance), so bus 101 sits
electrically very close to the REF bus: the measured dVm/d(susceptance) there is
~5.8e-5 pu/pu. Both devices are sized well past their PSS/E-typical rating (`number_of_steps`
/ `Y_increase` and `max_shunt_current`) so `|dy/dp|*(hi-lo)` clears `CONTROL_GAIN_FLOOR`
(src/definitions.jl) instead of being auto-frozen as insensitive
(src/discrete_control/control_continuation.jl:565-579), and each `voltage_setpoint`/
deadband is placed above the actual solved voltage (~1.0007 pu) so the continuation
actually drives the device toward its rail across several passes rather than finding it
already in-band on pass 1."""
function build_lcc_control_system(; p_set_mw::Union{Nothing, Float64} = nothing)
    raw = joinpath(TEST_DATA_DIR, "case5_2_lcc.raw")
    sys = make_system(PFP.PowerModelsData(raw); runchecks = false)
    bus101 = get_bus(sys, 101)
    add_component!(
        sys,
        SwitchedAdmittance(; name = "ctrl_shunt_101", available = true,
            bus = bus101, Y = 0.0 + 0.0im, initial_status = [0], number_of_steps = [8],
            Y_increase = [0.0 + 0.5im], admittance_limits = (min = 1.05, max = 1.08),
            control_mode = PSY.SwitchedAdmittanceControlMode.DISCRETE_VOLTAGE,
        ),
    )
    # FACTS construction mirrors build_ieee14_facts_system (test_discrete_control.jl:10-60)
    # with bus/name adjusted and `max_shunt_current`/`voltage_setpoint` widened/offset for the
    # same insensitive-bus reason as the shunt above (see docstring).
    add_component!(
        sys,
        FACTSControlDevice(;
            name = "ctrl_facts_101",
            available = true,
            bus = bus101,
            control_mode = PSY.FACTSOperationModes.NML,
            voltage_setpoint = 1.06,
            max_shunt_current = 1000.0,
            max_reactive_power = 9999.0,
            shunt_control_type = PSY.FACTSShuntControlType.STATCOM,
            regulated_bus_number = 0,
        ),
    )
    if p_set_mw !== nothing
        # `data.lcc.p_set` is seeded in `initialize_LCCParameters!` (src/lcc_utils.jl:861)
        # from `PSY.get_transfer_setpoint.(lccs)` in MW (divided by base power there to
        # reach pu), so the matching setter takes MW directly, not a pu value.
        for l in get_components(TwoTerminalLCCLine, sys)
            set_transfer_setpoint!(l, p_set_mw)
        end
    end
    return sys
end

@testset "LCC + control: equivalence with devices locked at converged settings" begin
    sys = build_lcc_control_system()
    pf = ACPolarPowerFlow(; control_discrete_devices = true)
    data = PowerFlowData(pf, sys)
    solve_power_flow!(data)
    @test all(data.converged)
    # The probe phase guarantees capture->solve->restore ran for each enrolled device.
    @test PowerFlows.get_control_inner_solve_count(data) >= 3

    results = PowerFlows.get_controlled_device_results(data)
    vm_ctrl = copy(data.bus_magnitude[:, 1])
    va_ctrl = copy(data.bus_angles[:, 1])
    lcc_taps = copy(data.lcc.rectifier.tap[:, 1])

    # Twin system: same devices fixed at the converged settings, control off.
    #   - shunt: `PSY.set_Y!(sa, Complex(d.g0, d.current))`, exactly what write_device_settings!
    #     (the continuation's save-back, src/solve_ac_power_flow.jl:130-181) writes -- this
    #     fixture's shunt is built with an all-zero `initial_status`, so `_is_binit_shunt`
    #     classifies it `psse_convention = true`
    #     (src/discrete_control/control_metadata.jl:147-149,327), the branch
    #     write_device_settings! takes; `d.g0 == real(get_Y(...)) == 0.0` here and `d.current`
    #     is exactly `results.final` for this row. A `SwitchedAdmittance`'s `Y` feeds the Ybus
    #     the same way whether or not `control_discrete_devices` is set, so this reproduces the
    #     controlled solve's physics exactly.
    #   - FACTS: write_device_settings! itself writes `PSY.set_reactive_power_required!(fd,
    #     delivered_q_mvar(d, |V|))`, and this test calls that same setter for reporting parity
    #     -- but it is a NO-OP for the actual physics: `contributes_reactive_power(::
    #     PSY.FACTSControlDevice) = false` (src/psi_utils.jl:61) makes `_get_injections!`
    #     (src/common.jl:20-40) skip every FACTSControlDevice unconditionally, so
    #     `reactive_power_required` is read only by reporting/export code (psse_export.jl),
    #     never by the AC solve. A `FACTSControlDevice`'s only physical footprint is the
    #     discrete-control continuation's live constant-Z withdrawal term
    #     (src/discrete_control/control_continuation.jl:297-310), which does not exist when
    #     `control_discrete_devices = false`. Confirmed empirically: locking only via
    #     `set_reactive_power_required!` left bus 101's twin voltage off by ~5e-4 -- the same
    #     order as the FACTS device's entire contribution (`dVm/db * b_final`), i.e. the FACTS
    #     term was silently absent, not merely approximated. The physically-equivalent lock is
    #     therefore a real shunt admittance of `results.final` (susceptance, not
    #     `delivered_q_mvar`) at the same bus, which -- like the `SwitchedAdmittance` above --
    #     feeds the Ybus the same way regardless of `control_discrete_devices`.
    sys2 = build_lcc_control_system()
    shunt_final = only(results.final[results.family .== "SwitchedAdmittance"])
    facts_final_b = only(results.final[results.family .== "FACTSControlDevice"])
    facts_q = only(results.delivered_q_mvar[results.family .== "FACTSControlDevice"])
    sa_orig = get_component(SwitchedAdmittance, sys, "ctrl_shunt_101")
    sa2 = get_component(SwitchedAdmittance, sys2, "ctrl_shunt_101")
    set_Y!(sa2, Complex(real(get_Y(sa_orig)), shunt_final))
    fd2 = get_component(FACTSControlDevice, sys2, "ctrl_facts_101")
    set_reactive_power_required!(fd2, facts_q)   # reporting parity only; see note above
    add_component!(
        sys2,
        FixedAdmittance(;
            name = "ctrl_facts_101_locked_shunt", available = true,
            bus = get_bus(sys2, 101), Y = Complex(0.0, facts_final_b),
        ),
    )

    pf2 = ACPolarPowerFlow(; control_discrete_devices = false)
    data2 = PowerFlowData(pf2, sys2)
    solve_power_flow!(data2)
    @test all(data2.converged)
    @test isapprox(data2.bus_magnitude[:, 1], vm_ctrl; atol = 1e-6)
    @test isapprox(data2.bus_angles[:, 1], va_ctrl; atol = 1e-6)
    @test isapprox(data2.lcc.rectifier.tap[:, 1], lcc_taps; atol = 1e-6)
end

@testset "LCC + control with 0 MW transfer (i_dc = 0 tap pinning) survives restores" begin
    sys = build_lcc_control_system(; p_set_mw = 0.0)
    pf = ACPolarPowerFlow(; control_discrete_devices = true)
    data = PowerFlowData(pf, sys)
    solve_power_flow!(data)
    @test all(data.converged)
    @test all(isapprox.(
        data.lcc.rectifier.tap[:, 1], data.lcc.rectifier.tap_setpoint; atol = 1e-8))
end
