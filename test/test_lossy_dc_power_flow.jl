@testset "Lossy DCLF: flat start produces same result as lossless" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    # Force flat-start: V=1, theta=0 on all buses — all computed losses should be zero.
    for bus in PSY.get_components(PSY.ACBus, sys)
        PSY.set_magnitude!(bus, 1.0)
        PSY.set_angle!(bus, 0.0)
    end
    lossless = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    lossy = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true, loss_approximation_as_injection = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    @test lossless["1"]["flow_results"].P_from_to ==
          lossy["1"]["flow_results"].P_from_to
    @test lossless["1"]["bus_results"].θ == lossy["1"]["bus_results"].θ
end

@testset "Lossy DCLF: default matches explicit false" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    default_result = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    explicit_false = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true, loss_approximation_as_injection = false),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    @test default_result["1"]["flow_results"].P_from_to ==
          explicit_false["1"]["flow_results"].P_from_to
end

@testset "Lossy DCLF: after AC solve, results differ from lossless" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    # Solve AC to populate a non-trivial voltage profile.
    solve_and_store_power_flow!(ACPowerFlow(; correct_bustypes = true), sys)

    lossless = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    lossy = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true, loss_approximation_as_injection = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    # Flows should differ (system has nonzero angles and V != 1.0).
    @test lossless["1"]["flow_results"].P_from_to !=
          lossy["1"]["flow_results"].P_from_to
    # Lossy model reports asymmetric flows: P_to_from != -P_from_to
    lossy_flows = lossy["1"]["flow_results"]
    @test lossy_flows.P_to_from != -lossy_flows.P_from_to
    # Loss conservation: P_from_to + P_to_from ≈ P_losses
    @test isapprox(
        lossy_flows.P_from_to .+ lossy_flows.P_to_from,
        lossy_flows.P_losses;
        atol = 1e-6,
    )
    # Lossless model still has symmetric flows.
    lossless_flows = lossless["1"]["flow_results"]
    @test lossless_flows.P_to_from == -lossless_flows.P_from_to
end

@testset "Lossy DCLF: initial_loss_injections field" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # Without loss_approximation_as_injection: field should be nothing.
    data_no_loss = PowerFlowData(DCPowerFlow(; correct_bustypes = true), sys)
    @test data_no_loss.initial_loss_injections === nothing

    # With loss_approximation_as_injection on flat-start: field should be all zeros.
    for bus in PSY.get_components(PSY.ACBus, sys)
        PSY.set_magnitude!(bus, 1.0)
        PSY.set_angle!(bus, 0.0)
    end
    data_flat = PowerFlowData(
        DCPowerFlow(; correct_bustypes = true, loss_approximation_as_injection = true), sys)
    @test data_flat.initial_loss_injections !== nothing
    @test all(data_flat.initial_loss_injections .== 0.0)

    # After AC solve: field should be non-zero.
    solve_and_store_power_flow!(ACPowerFlow(; correct_bustypes = true), sys)
    data_ac = PowerFlowData(
        DCPowerFlow(; correct_bustypes = true, loss_approximation_as_injection = true), sys)
    @test data_ac.initial_loss_injections !== nothing
    @test !all(data_ac.initial_loss_injections .== 0.0)
    # Loss injections should be non-positive (withdrawals).
    @test all(data_ac.initial_loss_injections .<= 0.0)
end

@testset "Lossy DCLF: hand-verified 3-bus system" begin
    sys = System(100.0)
    # Bus 1: REF, V=1.05, theta=0.0
    b1 = _add_simple_bus!(sys, 1, PSY.ACBusTypes.REF, 230, 1.05, 0.0)
    # Bus 2: PV, V=1.02, theta=-0.05
    b2 = _add_simple_bus!(sys, 2, PSY.ACBusTypes.PV, 230, 1.02, -0.05)
    # Bus 3: PQ, V=0.98, theta=-0.10
    b3 = _add_simple_bus!(sys, 3, PSY.ACBusTypes.PQ, 230, 0.98, -0.10)
    _add_simple_source!(sys, b1, 1.0, 0.0)
    _add_simple_thermal_standard!(sys, b2, 0.4, 0.0)
    _add_simple_load!(sys, b2, 0.3, 0.0)
    _add_simple_load!(sys, b3, 0.5, 0.0)
    _add_simple_line!(sys, b1, b2, 0.01, 0.1, 0.0)   # r=0.01, x=0.1
    _add_simple_line!(sys, b2, b3, 0.02, 0.15, 0.0)   # r=0.02, x=0.15
    _add_simple_line!(sys, b1, b3, 0.015, 0.12, 0.0)   # r=0.015, x=0.12

    data = PowerFlowData(
        DCPowerFlow(; loss_approximation_as_injection = true, correct_bustypes = true), sys)
    loss_inj = data.initial_loss_injections
    @test loss_inj !== nothing

    # Hand-compute expected losses using I²r where I = |Vi∠θi - Vj∠θj| / |Z|.
    # For lines (tap=1): |Z_eff|² = r² + x², |ΔV|² = Vi²+Vj²−2·Vi·Vj·cos(θi−θj).
    # P_loss = |ΔV|² · r / |Z_eff|²

    # Line 1→2: r=0.01, x=0.1, Vi=1.05, Vj=1.02, θi=0.0, θj=-0.05
    r12, x12 = 0.01, 0.1
    dV_sq_12 = 1.05^2 + 1.02^2 - 2 * 1.05 * 1.02 * cos(0.0 - (-0.05))
    P_loss_12 = dV_sq_12 * r12 / (r12^2 + x12^2)

    # Line 2→3: r=0.02, x=0.15, Vi=1.02, Vj=0.98, θi=-0.05, θj=-0.10
    r23, x23 = 0.02, 0.15
    dV_sq_23 = 1.02^2 + 0.98^2 - 2 * 1.02 * 0.98 * cos(-0.05 - (-0.10))
    P_loss_23 = dV_sq_23 * r23 / (r23^2 + x23^2)

    # Line 1→3: r=0.015, x=0.12, Vi=1.05, Vj=0.98, θi=0.0, θj=-0.10
    r13, x13 = 0.015, 0.12
    dV_sq_13 = 1.05^2 + 0.98^2 - 2 * 1.05 * 0.98 * cos(0.0 - (-0.10))
    P_loss_13 = dV_sq_13 * r13 / (r13^2 + x13^2)

    total_loss = P_loss_12 + P_loss_23 + P_loss_13
    # Total loss injections should equal negative of total losses.
    @test isapprox(sum(loss_inj), -total_loss, atol = 1e-10)

    # Sending ends: bus 1 (theta=0.0) sends to bus 2 (theta=-0.05) and bus 3 (theta=-0.10).
    # Bus 2 (theta=-0.05) sends to bus 3 (theta=-0.10).
    # So bus 1 receives losses from lines 1-2 and 1-3, bus 2 from line 2-3.
    bus_lookup = PF.get_bus_lookup(data)
    ix1 = bus_lookup[1]
    ix2 = bus_lookup[2]
    ix3 = bus_lookup[3]
    @test isapprox(loss_inj[ix1, 1], -(P_loss_12 + P_loss_13), atol = 1e-10)
    @test isapprox(loss_inj[ix2, 1], -P_loss_23, atol = 1e-10)
    @test loss_inj[ix3, 1] == 0.0
end

@testset "Lossy DCLF: multi-period applies same loss injection to all timesteps" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    solve_and_store_power_flow!(ACPowerFlow(; correct_bustypes = true), sys)
    data = PowerFlowData(
        DCPowerFlow(;
            loss_approximation_as_injection = true,
            time_steps = 3,
            correct_bustypes = true,
        ),
        sys,
    )
    @test data.initial_loss_injections !== nothing
    # All columns should be identical (same loss injection for each timestep).
    @test data.initial_loss_injections[:, 1] == data.initial_loss_injections[:, 2]
    @test data.initial_loss_injections[:, 1] == data.initial_loss_injections[:, 3]
end

@testset "Lossy DCLF: AC and PTDF data have nothing for initial_loss_injections" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    ac_data = PowerFlowData(ACPowerFlow(; correct_bustypes = true), sys)
    @test ac_data.initial_loss_injections === nothing
    ptdf_data = PowerFlowData(PTDFDCPowerFlow(; correct_bustypes = true), sys)
    @test ptdf_data.initial_loss_injections === nothing
end
