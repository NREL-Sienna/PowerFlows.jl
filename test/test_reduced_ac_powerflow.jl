ac_reduction_types = Dict{String, Vector{PNM.NetworkReduction}}(
    "default" => PNM.NetworkReduction[],
    "radial" => PNM.NetworkReduction[PNM.RadialReduction()],
    "degree 2" => PNM.NetworkReduction[PNM.DegreeTwoReduction()],
    "radial + degree 2" =>
        PNM.NetworkReduction[PNM.RadialReduction(), PNM.DegreeTwoReduction()],
)
@testset "AC power flow on 2k bus system: validate reduce-then-solve" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
    unreduced = PF.PowerFlowData(
        PF.ACPowerFlow(PF.TrustRegionACPowerFlow),
        sys;
        correct_bustypes = true,
    )
    PF.solve_powerflow!(unreduced)
    @assert all(unreduced.converged)
    pf = ACPowerFlow(PF.TrustRegionACPowerFlow)
    for (k, v) in ac_reduction_types
        size(v) == 0 && continue # no reduction at all.
        @testset "$k reduction" begin
            validate_reduced_powerflow(pf, sys, v, unreduced)
        end
    end
end

@testset "all reductions on psse_14_network_reduction_test_system" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    pf = ACPowerFlow(PF.TrustRegionACPowerFlow)

    for (k, v) in ac_reduction_types
        @testset "$k reduction" begin
            result = test_reduced_powerflow(pf, sys, v)
            @test all(result.converged)
        end
    end

    @testset "ward reduction" begin
        study_buses = [101, 114, 110, 111]
        result = test_reduced_powerflow(
            pf,
            sys,
            PNM.NetworkReduction[PNM.WardReduction(study_buses)],
        )
        @test all(result.converged) broken = true
    end
end
