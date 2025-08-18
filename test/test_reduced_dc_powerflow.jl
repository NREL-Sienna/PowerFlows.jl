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
        @testset "$dc_pf power flow" for dc_pf in dc_powerflows
            unreduced = PF.PowerFlowData(dc_pf, sys; correct_bustypes = true)
            PF.solve_powerflow!(unreduced)
            @assert all(unreduced.converged)
            validate_reduced_powerflow(dc_pf, sys, v, unreduced)
        end
    end
end
