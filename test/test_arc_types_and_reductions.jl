"""
Run one DC and one AC power flow on a system with a given set of network reductions.
"""
function test_all_power_flow_types(
    sys::System,
    network_reductions::Vector{PNM.NetworkReduction},
)
    # AC power flow has different syntax
    pf = ACPowerFlow{PF.TrustRegionACPowerFlow}()
    data = PF.PowerFlowData(
        pf,
        sys;
        network_reductions = deepcopy(network_reductions),
        correct_bustypes = true,
    )
    solve_power_flow!(data; pf = pf) # should run without errors.
    @test !isempty(data.arc_active_power_flow_from_to)
    solve_and_store_power_flow!(
        pf,
        sys;
        network_reductions = deepcopy(network_reductions),
        correct_bustypes = true,
    )
    pf = PF.DCPowerFlow()
    data = PF.PowerFlowData(pf, sys; network_reductions = deepcopy(network_reductions))
    solve_power_flow!(data) # should run without errors.
    @test !isempty(data.arc_active_power_flow_from_to)
end

@testset "radial reduction with 3WT/parallel lines" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[RadialReduction()]
    test_all_power_flow_types(sys, network_reductions)
end

@testset "degreee 2 reduction with 3WT/parallel lines" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[DegreeTwoReduction()]
    test_all_power_flow_types(sys, network_reductions)
end
