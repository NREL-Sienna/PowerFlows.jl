dc_powerflows = (
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
            PF.solve_powerflow!(unreduced)
            @assert all(unreduced.converged)
            validate_reduced_powerflow(dc_pf, sys, v, unreduced)
            results = PF.solve_powerflow(
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
