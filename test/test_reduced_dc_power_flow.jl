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
            dc_pf = pf_type(; correct_bustypes = true)
            unreduced = PF.PowerFlowData(dc_pf, sys)
            PF.solve_power_flow!(unreduced)
            @assert all(unreduced.converged)
            dc_pf_reduced =
                pf_type(; network_reductions = deepcopy(v), correct_bustypes = true)
            validate_reduced_power_flow(dc_pf_reduced, sys, v, unreduced)
            results = PF.solve_power_flow(dc_pf_reduced, sys, PF.FlowReporting.ARC_FLOWS)
            # no write-results-to-system solve_and_store_power_flow! for DC: should we add one?
        end
    end
end
