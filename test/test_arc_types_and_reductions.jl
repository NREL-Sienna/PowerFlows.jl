function test_all_powerflow_types(
    sys::System,
    network_reductions::Vector{PNM.NetworkReduction},
)
    # AC powerflow has different syntax
    pf = ACPowerFlow()
    data = PF.PowerFlowData(
        pf,
        sys;
        network_reductions = deepcopy(network_reductions),
        correct_bustypes = true,
    )
    solve_powerflow!(data; pf = pf) # should run without errors.
    @test !isempty(data.arc_activepower_flow_from_to)
    solve_powerflow!(
        pf,
        sys;
        network_reductions = deepcopy(network_reductions),
        correct_bustypes = true,
    )
    pf = PF.DCPowerFlow()
    data = PF.PowerFlowData(pf, sys; network_reductions = deepcopy(network_reductions))
    solve_powerflow!(data) # should run without errors.
    @test !isempty(data.arc_activepower_flow_from_to)
end

@testset "radial reduction with 3WT/parallel lines" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[RadialReduction()]
    test_all_powerflow_types(sys, network_reductions)
end

@testset "degreee 2 reduction with 3WT/parallel lines" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[DegreeTwoReduction()]
    test_all_powerflow_types(sys, network_reductions)
end
