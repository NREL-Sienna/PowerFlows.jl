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

@testset "system + powerflow solver calls" begin
    for (k, v) in ac_reduction_types
        @testset "$k reduction" begin
            # systems in PSB with 3WT's:
            #=
            (PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
            (PSSEParsingTestSystems, "psse_14_tap_correction_test_system")
            (PSSEParsingTestSystems, "psse_14_zero_impedance_branch_test_system")
            (PSSEParsingTestSystems, "psse_4_zero_impedance_3wt_test_system")
            (PSSEParsingTestSystems, "psse_ybus_14_test_system")
            (PSSEParsingTestSystems, "pti_case10_voltage_winding_correction_sys")
            (PSSEParsingTestSystems, "pti_case8_voltage_winding_correction_sys")
            (PSSEParsingTestSystems, "pti_frankenstein_20_sys")
            (PSSEParsingTestSystems, "pti_frankenstein_70_sys")
            (PSSEParsingTestSystems, "pti_modified_case14_sys")
            (PSSEParsingTestSystems, "pti_three_winding_mag_test_sys")
            (PSSEParsingTestSystems, "pti_three_winding_test_2_sys")
            (PSSEParsingTestSystems, "pti_three_winding_test_sys")
            =#
            sys = build_system(
                PSSEParsingTestSystems,
                "psse_14_network_reduction_test_system",
            )
            pf = ACPowerFlow(PF.TrustRegionACPowerFlow)
            results = solve_powerflow(
                pf,
                sys;
                correct_bustypes = true,
                network_reductions = deepcopy(v),
            )
            @assert !isempty(PSY.get_components(PSY.Transformer3W, sys))
            test_trf = first(collect(PSY.get_components(PSY.Transformer3W, sys)))
            test_trf_name = PSY.get_name(test_trf)
            if !ismissing(results)
                arc_flows = results["flow_results"]
                trf_arc_flows = zeros(ComplexF32, 3)
                for i in 1:3
                    adj = ("primary", "secondary", "tertiary")[i]
                    ix = arc_flows[!, "line_name"] .== "$(test_trf_name)-$adj"
                    trf_arc_flows[i] =
                        sum(arc_flows[ix, "P_from_to"]) +
                        im * sum(arc_flows[ix, "Q_from_to"])
                end
                @test solve_powerflow!(
                    pf,
                    sys;
                    correct_bustypes = true,
                    network_reductions = deepcopy(v),
                )
                base_power = PSY.get_base_power(sys)
                # check that transformer bus-to-star entries are there.
                @test isapprox(
                    PSY.get_active_power_flow_primary(test_trf),
                    real(trf_arc_flows[1]) / base_power;
                    atol = 1e-5,
                )
                @test isapprox(
                    PSY.get_reactive_power_flow_primary(test_trf),
                    imag(trf_arc_flows[1]) / base_power;
                    atol = 1e-5,
                )
                @test isapprox(
                    PSY.get_active_power_flow_secondary(test_trf),
                    real(trf_arc_flows[2]) / base_power;
                    atol = 1e-5,
                )
                @test isapprox(
                    PSY.get_reactive_power_flow_secondary(test_trf),
                    imag(trf_arc_flows[2]) / base_power;
                    atol = 1e-5,
                )
                @test isapprox(
                    PSY.get_active_power_flow_tertiary(test_trf),
                    real(trf_arc_flows[3]) / base_power;
                    atol = 1e-5,
                )
                @test isapprox(
                    PSY.get_reactive_power_flow_tertiary(test_trf),
                    imag(trf_arc_flows[3]) / base_power;
                    atol = 1e-5,
                )
            else
                @error "Failed to converge: cannot test AC post-processing with $k reduction"
            end
        end
    end
end
