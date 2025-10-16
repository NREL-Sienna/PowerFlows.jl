const UNSUPPORTED =
    Set{Tuple{Type{<:PNM.NetworkReduction}, Type{<:PF.PowerFlowEvaluationModel}}}(
        [
        (PNM.WardReduction, PF.ACPowerFlow),
        (PNM.WardReduction, PF.DCPowerFlow),
        (PNM.WardReduction, PF.PTDFDCPowerFlow),
        (PNM.WardReduction, PF.vPTDFDCPowerFlow),
        (PNM.RadialReduction, PF.ACPowerFlow),
    ],
    )

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
        if any((typeof(nr), typeof(pf)) in UNSUPPORTED for nr in v)
            @warn "Skipping unsupported combination: $(typeof(nr)), $(typeof(pf))"
            continue
        end
        @testset "$k reduction" begin
            validate_reduced_powerflow(pf, sys, v, unreduced)
        end
    end
end

@testset "all reductions on psse_14_network_reduction_test_system" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    pf = ACPowerFlow(PF.TrustRegionACPowerFlow)

    for (k, v) in ac_reduction_types
        if any((typeof(nr), typeof(pf)) in UNSUPPORTED for nr in v)
            @warn "Skipping unsupported combination: $(typeof(nr)), $(typeof(pf))"
            continue
        end
        @testset "$k reduction" begin
            result = test_reduced_powerflow(pf, sys, v)
            @test all(result.converged)
        end
    end

    # not yet implemented.
    #=
    @testset "ward reduction" begin
        study_buses = [101, 114, 110, 111]
        result = test_reduced_powerflow(
            pf,
            sys,
            PNM.NetworkReduction[PNM.WardReduction(study_buses)],
        )
        @test all(result.converged) broken = true
    end
    =#
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

function compare_voltages(
    unreduced::PF.PowerFlowData,
    sys::PSY.System,
    temp_bus_map::Dict{Int, String},
    bus_no::Int,
)
    bus_lookup = PF.get_bus_lookup(unreduced)
    unreduced_Vm = unreduced.bus_magnitude[bus_lookup[bus_no], 1]
    unreduced_Va = unreduced.bus_angles[bus_lookup[bus_no], 1]
    bus_name = temp_bus_map[bus_no]
    bus = PSY.get_component(PSY.ACBus, sys, bus_name)
    reduced_Vm = PSY.get_magnitude(bus)
    reduced_Va = PSY.get_angle(bus)
    @test isapprox(unreduced_Vm, reduced_Vm; atol = 1e-6)
    @test isapprox(unreduced_Va, reduced_Va; atol = 1e-6)
end

function compare_power_flows(
    unreduced::PF.PowerFlowData,
    sys::PSY.System,
    branch::PSY.Branch,
)
    name = PSY.get_name(branch)
    arc_lookup = PF.get_arc_lookup(unreduced)
    arc_ix = arc_lookup[PNM.get_arc_tuple(branch)]
    unreduced_active_flow = unreduced.arc_activepower_flow_from_to[arc_ix, 1]
    unreduced_reactive_flow = unreduced.arc_reactivepower_flow_from_to[arc_ix, 1]
    reduced_active_flow =
        PSY.get_active_power_flow(PSY.get_component(PSY.Branch, sys, name))
    reduced_reactive_flow =
        PSY.get_reactive_power_flow(PSY.get_component(PSY.Branch, sys, name))
    @test isapprox(unreduced_active_flow, reduced_active_flow; atol = 1e-3)
    @test isapprox(unreduced_reactive_flow, reduced_reactive_flow; atol = 1e-3)
end

@testset "parallel branches: recovering flows" begin
    # we can't run an unreduced power flow and compare: parallel branches are always reduced.
    # instead, do a sanity check: flows over the parallel branches should add up to
    # the net flow across the equivalent branch.

    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
    unreduced = PF.PowerFlowData(
        PF.ACPowerFlow(PF.TrustRegionACPowerFlow),
        sys;
        correct_bustypes = true,
    )
    pf = PF.ACPowerFlow{PF.TrustRegionACPowerFlow}(; skip_redistribution = true)
    PF.solve_powerflow!(unreduced; pf = pf)
    PF.solve_powerflow!(
        pf,
        sys;
        correct_bustypes = true,
        network_reductions = PNM.NetworkReduction[PNM.DegreeTwoReduction()],
    )
    temp_ybus =
        PNM.Ybus(sys; network_reductions = PNM.NetworkReduction[PNM.DegreeTwoReduction()])
    nrd = PNM.get_network_reduction_data(temp_ybus)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )
    parallel_br_map = PNM.get_parallel_branch_map(nrd)
    arc_lookup = PF.get_arc_lookup(unreduced)
    for (equiv_arc, branches) in parallel_br_map
        equiv_arc_ix = arc_lookup[equiv_arc]
        net_flow_from_to =
            unreduced.arc_activepower_flow_from_to[equiv_arc_ix, 1] +
            im * unreduced.arc_reactivepower_flow_from_to[equiv_arc_ix, 1]
        net_flow_to_from =
            unreduced.arc_activepower_flow_to_from[equiv_arc_ix, 1] +
            im * unreduced.arc_reactivepower_flow_to_from[equiv_arc_ix, 1]
        total_flow = zero(ComplexF32)
        expected_from_bus = equiv_arc[1]
        (from_bus_no, to_bus_no) = PNM.get_arc_tuple(first(branches))
        @assert equiv_arc == (from_bus_no, to_bus_no) ||
                equiv_arc == (to_bus_no, from_bus_no)
        reversed = from_bus_no != expected_from_bus
        for br in branches
            @assert PNM.get_arc_tuple(br) == (from_bus_no, to_bus_no)
            total_flow +=
                PSY.get_active_power_flow(br) +
                im * PSY.get_reactive_power_flow(br)
        end
        if reversed
            @test isapprox(net_from_to_from, total_flow; atol = 1e-3)
        else
            @test isapprox(net_flow_from_to, total_flow; atol = 1e-3)
        end
    end
end

@testset "degree 2 reduction: recovering flows/voltages" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
    unreduced = PF.PowerFlowData(
        PF.ACPowerFlow(PF.TrustRegionACPowerFlow),
        sys;
        correct_bustypes = true,
    )
    pf = PF.ACPowerFlow{PF.TrustRegionACPowerFlow}(; skip_redistribution = true)
    PF.solve_powerflow!(unreduced; pf = pf)
    PF.solve_powerflow!(
        pf,
        sys;
        correct_bustypes = true,
        network_reductions = PNM.NetworkReduction[PNM.DegreeTwoReduction()],
    )
    temp_ybus =
        PNM.Ybus(sys; network_reductions = PNM.NetworkReduction[PNM.DegreeTwoReduction()])
    nrd = PNM.get_network_reduction_data(temp_ybus)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )
    bus_lookup = PF.get_bus_lookup(unreduced)
    for (equivalent_arc, segments) in PNM.get_series_branch_map(nrd)
        for segment in segments
            (from_bus_no, to_bus_no) = PNM.get_arc_tuple(segment)
            compare_voltages(unreduced, sys, temp_bus_map, from_bus_no)
            compare_voltages(unreduced, sys, temp_bus_map, to_bus_no)
            if !(from_bus_no in equivalent_arc)
                @assert unreduced.bus_type[bus_lookup[from_bus_no]] == PSY.ACBusTypes.PQ
            end
            if !(to_bus_no in equivalent_arc)
                @assert unreduced.bus_type[bus_lookup[to_bus_no]] == PSY.ACBusTypes.PQ
            end
        end
    end
    for (equivalent_arc, segments) in PNM.get_series_branch_map(nrd)
        for segment in segments
            # skip parallel branches
            if segment isa PSY.ACTransmission
                compare_power_flows(unreduced, sys, segment)
            end
        end
    end
end
