function test_reduced_ac_powerflow(nrs::Vector{PNM.NetworkReduction})
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    pf = PF.ACPowerFlow(PF.NewtonRaphsonACPowerFlow)
    data = PF.PowerFlowData(pf, sys; time_steps = 1, network_reductions = nrs)
    PF.solve_powerflow!(data; pf = pf)
    return data.converged
end

@testset "14 bus with default reductions" begin
    isconverged = test_reduced_ac_powerflow(PNM.NetworkReduction[])
    @test all(isconverged)
end

@testset "14 bus with radial reduction" begin
    isconverged = test_reduced_ac_powerflow(PNM.NetworkReduction[PNM.RadialReduction()])
    @test all(isconverged)
end

@testset "14 bus with degree 2 reduction" begin
    isconverged = test_reduced_ac_powerflow(PNM.NetworkReduction[PNM.DegreeTwoReduction()])
    @test all(isconverged)
end

@testset "14 bus with ward reduction" begin
    study_buses = [101, 114, 110, 111]
    isconverged =
        test_reduced_ac_powerflow(PNM.NetworkReduction[PNM.WardReduction(study_buses)])
    @test all(isconverged)
end
