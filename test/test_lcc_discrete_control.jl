# Tests for LCC coexistence with the discrete-control continuation.
#
# The classification below is the contract the continuation's checkpoint implements: every
# per-time-step Matrix field of the LCC state belongs to exactly one bucket, so a new field
# added without extending the checkpoint fails a test rather than breaking rollback silently.
# `i_dc` is CHECKPOINTED because `initialize_LCCParameters!` derives it from `p_set` at
# construction and no solve path rewrites it.

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

    # Driven off the classification rather than naming columns literally, so a newly
    # classified field is exercised here automatically.
    checkpointed_cols = vcat(
        [
            (conv, f)
            for conv in (data.lcc.rectifier, data.lcc.inverter)
            for f in LCC_CHECKPOINTED_FIELDS[PF.LCCConverterParameters]
        ],
        [(data.lcc, f) for f in LCC_CHECKPOINTED_FIELDS[PF.LCCParameters]],
    )

    # Establish a consistent derived state, then snapshot.
    PowerFlows._update_ybus_lcc!(data, ts)
    ref_cols = [copy(getfield(obj, f)[:, ts]) for (obj, f) in checkpointed_cols]
    ref_admittances = copy(data.lcc.branch_admittances)
    snap = PowerFlows._snapshot_state(data, ts)
    # 8 non-LCC columns plus one per checkpointed column: catches `_lcc_state_cols` drifting
    # from the classification.
    @test length(snap) == 8 + length(checkpointed_cols)

    # Simulate a diverged trial. Garbage is distinct per column so a column swap fails too.
    for (i, (obj, f)) in enumerate(checkpointed_cols)
        getfield(obj, f)[:, ts] .= 100.0 + i
    end
    data.bus_magnitude[:, ts] .= 0.5
    PowerFlows._update_ybus_lcc!(data, ts)  # caches now reflect the garbage

    PowerFlows._restore_state!(data, ts, snap)

    for ((obj, f), ref) in zip(checkpointed_cols, ref_cols)
        @test getfield(obj, f)[:, ts] == ref
    end
    # Caches recomputed at the restored state, not left holding the garbage-derived values.
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
switched shunt and a shunt FACTS device at bus 101 (PQ, 230 kV, largest load). The fixture
has no transformers, so tap devices stay covered by the existing discrete-control suites.
`p_set_mw` overrides both LCC transfer setpoints (0.0 exercises the i_dc = 0 tap-pinning
branch).

Every branch in the fixture carries x = 1e-4 pu, leaving bus 101 electrically bolted to the
REF bus (dVm/d(susceptance) there is ~5.8e-5 pu/pu). Device ratings and setpoints are
therefore sized far past anything realistic: the ratings clear `CONTROL_GAIN_FLOOR` so the
devices enroll instead of being frozen as insensitive, and the setpoints sit above the
solved voltage (~1.0007 pu) so the continuation drives them over several passes."""
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
        # `initialize_LCCParameters!` seeds `p_set` from this setter's value in MW.
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

    # Twin system: same devices fixed at the converged settings, control off. The shunt locks
    # through its own `Y`, as `write_device_settings!` saves back. The FACTS device needs a
    # `FixedAdmittance` instead: `contributes_reactive_power(::FACTSControlDevice)` is false,
    # so `reactive_power_required` reaches reporting only, and the device's physical footprint
    # is a constant-Z withdrawal that exists solely inside the continuation. Locking via
    # `set_reactive_power_required!` alone leaves the twin's voltage off by the device's whole
    # contribution (~5e-4 pu here).
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
    @test all(
        isapprox.(
            data.lcc.rectifier.tap[:, 1], data.lcc.rectifier.tap_setpoint; atol = 1e-8),
    )
end

"""Per-step reactive-withdrawal scale for bus 101, mirroring `_shunt_step_q_scale`
(test/test_utils/common.jl). Scaling is narrowed to bus 101, the bus hosting the enrolled
devices, so the fixture's other loads stay fixed."""
_lcc_bus101_q_scale(t::Int) = 0.8 + 0.2 * t

"""Scale bus 101's reactive withdrawal across time-step columns `1:n` of a multiperiod
`data`, taking column 1 as the base."""
function _set_lcc_bus101_loads!(data, n::Int)
    bix = PF.get_bus_lookup(data)[101]
    base_q = data.bus_reactive_power_withdrawals[bix, 1]
    for t in 1:n
        data.bus_reactive_power_withdrawals[bix, t] = base_q * _lcc_bus101_q_scale(t)
    end
    return
end

"""Set bus 101's reactive withdrawal in a single-period `data` to the value
`_set_lcc_bus101_loads!` assigns to multiperiod column `t`, for building one step's twin."""
function _set_lcc_bus101_load_at_step!(data, t::Int)
    bix = PF.get_bus_lookup(data)[101]
    data.bus_reactive_power_withdrawals[bix, 1] *= _lcc_bus101_q_scale(t)
    return
end

@testset "LCC + control multiperiod: per-step equivalence" begin
    time_steps = 3
    sys = build_lcc_control_system()
    pf = ACPolarPowerFlow(; control_discrete_devices = true, time_steps = time_steps)
    data = PowerFlowData(pf, sys)
    _set_lcc_bus101_loads!(data, time_steps)
    solve_power_flow!(data)
    @test all(data.converged)
    # The per-step loads must move the solution clear of the 1e-8 tolerance the isolation
    # test below uses, or that test could pass with nothing left to tell the steps apart.
    @test maximum(abs.(data.bus_magnitude[:, 1] .- data.bus_magnitude[:, 3])) > 1e-7
    for ts in 1:time_steps
        # Twin lock as in the single-period equivalence test, per step.
        results = PowerFlows.get_controlled_device_results(data)
        step_results = results[results.time_step .== ts, :]
        shunt_final = only(step_results.final[step_results.family .== "SwitchedAdmittance"])
        facts_final_b =
            only(step_results.final[step_results.family .== "FACTSControlDevice"])
        facts_q =
            only(
                step_results.delivered_q_mvar[step_results.family .== "FACTSControlDevice"],
            )

        sys2 = build_lcc_control_system()
        sa_orig = get_component(SwitchedAdmittance, sys, "ctrl_shunt_101")
        sa2 = get_component(SwitchedAdmittance, sys2, "ctrl_shunt_101")
        set_Y!(sa2, Complex(real(get_Y(sa_orig)), shunt_final))
        fd2 = get_component(FACTSControlDevice, sys2, "ctrl_facts_101")
        set_reactive_power_required!(fd2, facts_q)   # reporting parity only, not physics
        add_component!(
            sys2,
            FixedAdmittance(;
                name = "ctrl_facts_101_locked_shunt", available = true,
                bus = get_bus(sys2, 101), Y = Complex(0.0, facts_final_b),
            ),
        )

        pf2 = ACPolarPowerFlow(; control_discrete_devices = false)
        data2 = PowerFlowData(pf2, sys2)
        _set_lcc_bus101_load_at_step!(data2, ts)
        solve_power_flow!(data2)
        @test all(data2.converged)
        @test isapprox(data2.bus_magnitude[:, 1], data.bus_magnitude[:, ts]; atol = 1e-6)
        @test isapprox(
            data2.lcc.rectifier.tap[:, 1], data.lcc.rectifier.tap[:, ts]; atol = 1e-6)
    end
end

@testset "LCC + control multiperiod: probe restores do not leak across steps" begin
    time_steps = 3
    sys = build_lcc_control_system()
    pf = ACPolarPowerFlow(; control_discrete_devices = true, time_steps = time_steps)
    data = PowerFlowData(pf, sys)
    _set_lcc_bus101_loads!(data, time_steps)
    solve_power_flow!(data)
    @test all(data.converged)
    # Independent single-period solves of steps 1 and 3 must match the multiperiod run's,
    # proving step 2's probes and rollbacks left no trace on its neighbours.
    for ts in (1, 3)
        sys_s = build_lcc_control_system()
        pf_s = ACPolarPowerFlow(; control_discrete_devices = true)
        data_s = PowerFlowData(pf_s, sys_s)
        _set_lcc_bus101_load_at_step!(data_s, ts)
        solve_power_flow!(data_s)
        @test all(data_s.converged)
        @test isapprox(data_s.bus_magnitude[:, 1], data.bus_magnitude[:, ts]; atol = 1e-8)
        @test isapprox(data_s.lcc.i_dc[:, 1], data.lcc.i_dc[:, ts]; atol = 1e-8)
    end
end

@testset "LCC + control multiperiod: single-step i_dc = 0 tap pinning" begin
    time_steps = 3
    # Unmodified reference run: identical construction/loads to the equivalence test above.
    sys_ref = build_lcc_control_system()
    pf = ACPolarPowerFlow(; control_discrete_devices = true, time_steps = time_steps)
    data_ref = PowerFlowData(pf, sys_ref)
    _set_lcc_bus101_loads!(data_ref, time_steps)
    solve_power_flow!(data_ref)
    @test all(data_ref.converged)

    # Same construction, but force LCC transfer to 0 MW at step 2 only.
    sys = build_lcc_control_system()
    data = PowerFlowData(pf, sys)
    _set_lcc_bus101_loads!(data, time_steps)
    # `i_dc` is derived from `p_set` at construction and never recomputed during a solve, so
    # overriding `p_set` afterwards requires zeroing `i_dc` too.
    data.lcc.p_set[:, 2] .= 0.0
    data.lcc.i_dc[:, 2] .= 0.0
    solve_power_flow!(data)
    @test all(data.converged)

    @test all(
        isapprox.(
            data.lcc.rectifier.tap[:, 2], data.lcc.rectifier.tap_setpoint; atol = 1e-8),
    )
    for ts in (1, 3)
        @test isapprox(
            data.bus_magnitude[:, ts],
            data_ref.bus_magnitude[:, ts];
            atol = 1e-6,
        )
        @test isapprox(
            data.lcc.rectifier.tap[:, ts], data_ref.lcc.rectifier.tap[:, ts]; atol = 1e-6)
    end
end

@testset "LCC + control: non-polar formulation (ACRectangularPowerFlow)" begin
    # The linear sensitivity path is polar-Jacobian-only, so a non-polar formulation drives
    # the continuation through the FD `_plant_sign` probe fallback instead.
    sys_polar = build_lcc_control_system()
    pf_polar = ACPolarPowerFlow(; control_discrete_devices = true)
    data_polar = PowerFlowData(pf_polar, sys_polar)
    solve_power_flow!(data_polar)
    @test all(data_polar.converged)

    sys_rect = build_lcc_control_system()
    pf_rect = ACRectangularPowerFlow(; control_discrete_devices = true)
    data_rect = PowerFlowData(pf_rect, sys_rect)
    solve_power_flow!(data_rect)
    @test all(data_rect.converged)

    # Rectangular-CI carries an e/f state internally but still populates `bus_magnitude`, so
    # the physical quantities compare directly.
    @test isapprox(
        data_rect.bus_magnitude[:, 1],
        data_polar.bus_magnitude[:, 1];
        atol = 1e-6,
    )
    @test isapprox(
        data_rect.lcc.rectifier.tap[:, 1], data_polar.lcc.rectifier.tap[:, 1]; atol = 1e-6)
    @test isapprox(
        data_rect.lcc.inverter.tap[:, 1], data_polar.lcc.inverter.tap[:, 1]; atol = 1e-6)
    @test isapprox(
        data_rect.lcc.rectifier.thyristor_angle[:, 1],
        data_polar.lcc.rectifier.thyristor_angle[:, 1];
        atol = 1e-6,
    )
    @test isapprox(data_rect.lcc.i_dc[:, 1], data_polar.lcc.i_dc[:, 1]; atol = 1e-6)
end
