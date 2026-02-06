@testset "DC power flow arc/branch flow reporting" begin
    sys = build_system(PSITestSystems, "case10_radial_series_reductions")
    pf_data_full = PF.PowerFlowData(
        DCPowerFlow(),
        sys,
    )
    pf_data_reduced = PF.PowerFlowData(
        DCPowerFlow(; network_reductions = NetworkReduction[DegreeTwoReduction()]),
        sys,
    )
    res_full_branch =
        solve_power_flow(pf_data_full, sys, PF.FlowReporting.BRANCH_FLOWS)["1"]["flow_results"]
    res_full_arc =
        solve_power_flow(pf_data_full, sys, PF.FlowReporting.ARC_FLOWS)["1"]["flow_results"]
    res_reduced_branch =
        solve_power_flow(pf_data_reduced, sys, PF.FlowReporting.BRANCH_FLOWS)["1"]["flow_results"]
    res_reduced_arc =
        solve_power_flow(pf_data_reduced, sys, PF.FlowReporting.ARC_FLOWS)["1"]["flow_results"]
    n_arcs = length(get_components(Arc, sys))
    @test size(res_full_arc)[1] == n_arcs

    two_arc_branches = length(
        get_available_components(
            x -> !isa(x, ThreeWindingTransformer),
            ACTransmission,
            sys,
        ),
    )
    three_arc_branches =
        length(get_components(get_available_primary, ThreeWindingTransformer, sys)) +
        length(get_components(get_available_secondary, ThreeWindingTransformer, sys)) +
        length(get_components(get_available_tertiary, ThreeWindingTransformer, sys))
    n_branches = two_arc_branches + three_arc_branches
    @test size(res_full_branch)[1] == n_branches

    @test size(res_reduced_arc)[1] < n_arcs
    @test size(res_reduced_branch)[1] == n_branches

    # This is the key test; we can recover exactly the same flow results dataframe after degree two reduction is used in the power flow calculation
    for (full_row, reduced_row) in
        zip(eachrow(res_full_branch), eachrow(res_reduced_branch))
        for name in names(full_row)
            if name == "flow_name"
                @test full_row[name] == reduced_row[name]
            else
                @test isapprox(full_row[name], reduced_row[name])
            end
        end
    end
end
