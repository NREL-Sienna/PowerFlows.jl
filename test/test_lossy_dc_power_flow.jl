@testset "Lossy DCLF: default matches explicit lossy_flows=false" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    default_result = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    explicit_false = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true, lossy_flows = false),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    @test default_result["1"]["flow_results"].P_from_to ==
          explicit_false["1"]["flow_results"].P_from_to
end

@testset "Lossy DCLF: lossy_flows=true gives different and asymmetric flows" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    lossless = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    lossy = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true, lossy_flows = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    # Flows should differ (c_sys14 has nonzero resistance).
    @test lossless["1"]["flow_results"].P_from_to != lossy["1"]["flow_results"].P_from_to
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

@testset "Lossy DCLF: arc_lossy_admittance fields" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # Without lossy_flows: arc admittance fields should be nothing.
    data_no_loss = PowerFlowData(DCPowerFlow(; correct_bustypes = true), sys)
    @test data_no_loss.arc_lossy_admittance_from_to === nothing
    @test data_no_loss.arc_lossy_admittance_to_from === nothing

    # With lossy_flows: arc admittance fields should be set.
    data_lossy =
        PowerFlowData(DCPowerFlow(; correct_bustypes = true, lossy_flows = true), sys)
    @test data_lossy.arc_lossy_admittance_from_to !== nothing
    @test data_lossy.arc_lossy_admittance_to_from !== nothing
end

@testset "Lossy DCLF: multi-period" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    data = PowerFlowData(
        DCPowerFlow(;
            lossy_flows = true,
            time_steps = 3,
            correct_bustypes = true,
        ),
        sys,
    )
    # Replicate time step 1 data (injections/withdrawals) to all time steps.
    for col in 2:3
        data.bus_active_power_injections[:, col] .=
            data.bus_active_power_injections[:, 1]
        data.bus_active_power_withdrawals[:, col] .=
            data.bus_active_power_withdrawals[:, 1]
    end
    solve_power_flow!(data)
    # All time steps should produce identical results (same system state).
    @test data.arc_active_power_flow_from_to[:, 1] ==
          data.arc_active_power_flow_from_to[:, 2]
    @test data.arc_active_power_flow_from_to[:, 1] ==
          data.arc_active_power_flow_from_to[:, 3]
end

@testset "Lossy DCLF: AC and PTDF data have no lossy admittances" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    ac_data = PowerFlowData(ACPowerFlow(; correct_bustypes = true), sys)
    @test ac_data.arc_lossy_admittance_from_to === nothing
    ptdf_data = PowerFlowData(PTDFDCPowerFlow(; correct_bustypes = true), sys)
    @test ptdf_data.arc_lossy_admittance_from_to === nothing
end
