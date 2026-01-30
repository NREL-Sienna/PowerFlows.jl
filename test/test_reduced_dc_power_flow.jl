dc_power_flows = (
    PF.DCPowerFlow(),
    PF.PTDFDCPowerFlow(),
    PF.vPTDFDCPowerFlow(),
)
dc_reduction_types = Dict{String, Vector{PNM.NetworkReduction}}(
    "radial" => PNM.NetworkReduction[PNM.RadialReduction()],
    "degree 2" => PNM.NetworkReduction[PNM.DegreeTwoReduction()],
    "radial + degree 2" =>
        PNM.NetworkReduction[PNM.RadialReduction(), PNM.DegreeTwoReduction()],
)
@testset "DC power flow on 2k bus system: validate reduce-then-solve" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
    @testset "$k reduction" for (k, v) in dc_reduction_types
        @testset "$pf_type power flow" for pf_type in
                                           InteractiveUtils.subtypes(PF.AbstractDCPowerFlow)
            dc_pf = pf_type()
            unreduced = PF.PowerFlowData(dc_pf, sys; correct_bustypes = true)
            PF.solve_power_flow!(unreduced)
            @assert all(unreduced.converged)
            validate_reduced_power_flow(dc_pf, sys, v, unreduced)
            results = PF.solve_power_flow(
                dc_pf,
                sys;
                network_reductions = deepcopy(v),
                correct_bustypes = true,
            )
            # no write-results-to-system solve_and_store_power_flow! for DC: should we add one?
            # PF.solve_and_store_power_flow!(dc_pf, sys; network_reductions = deepcopy(v), correct_bustypes = true)
        end
    end
end

# Ward reduction tested separately from Radial/Degree Two reductions because power flow results for retained buses is not expected to match exactly.
@testset "RTS: Ward Reduction and DC Power Flow" begin
    sys = build_system(PSISystems, "RTS_GMLC_DA_sys")
    bus_numbers = get_number.(get_components(ACBus, sys))
    study_buses = filter!(x -> digits(x)[end] == 1, bus_numbers)  #study buses are from area 1
    n_study_buses = length(study_buses)
    @testset "ward reduction $pf_type power flow" for pf_type in
                                                      InteractiveUtils.subtypes(
        PF.AbstractDCPowerFlow,
    )
        dc_pf = pf_type()
        unreduced = PF.PowerFlowData(dc_pf, sys; correct_bustypes = true)
        PF.solve_power_flow!(unreduced)
        @assert all(unreduced.converged)
        results = PF.solve_power_flow(
            dc_pf,
            sys;
            network_reductions = NetworkReduction[WardReduction(study_buses)],
            correct_bustypes = true,
        )
        @test nrow(results["1"]["bus_results"]) == n_study_buses
    end
end
