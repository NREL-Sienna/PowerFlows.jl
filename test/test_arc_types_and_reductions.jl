function types_in_series_reduction(nrd::PNM.NetworkReductionData)
    types = Set{DataType}()
    for segments in values(PNM.get_series_branch_map(nrd))
        for comp in segments
            push!(types, typeof(comp))
        end
    end
    return types
end

function find_parallel_arc(sys::System)
    arcs_seen = Set{Tuple{Int, Int}}()
    for br in PSY.get_components(PSY.ACBranch, sys)
        arc = PNM.get_arc_tuple(br)
        if arc in arcs_seen
            return arc
        else
            push!(arcs_seen, arc)
        end
    end
    error("No parallel arcs found")
    return (-1, -1)
end

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
    @test data isa Any
    solve_powerflow!(
        pf,
        sys;
        network_reductions = deepcopy(network_reductions),
        correct_bustypes = true,
    )
    for pf in [PF.DCPowerFlow(), PF.PTDFDCPowerFlow(), PF.vPTDFDCPowerFlow()]
        data = PF.PowerFlowData(pf, sys; network_reductions = deepcopy(network_reductions))
        solve_powerflow!(data) # should run without errors.
        @test data isa Any
    end
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
