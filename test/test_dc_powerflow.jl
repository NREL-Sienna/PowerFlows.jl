function flows_from_dataframe(powerflow_results,
    branch_lookup::Dict{Tuple{Int, Int}, Int},
    direction::Symbol = :P_from_to,
)
    flow_results_df = powerflow_results["1"]["flow_results"]
    flows = fill(NaN, length(branch_lookup))
    for row in eachrow(flow_results_df)
        flows[branch_lookup[(row.bus_from, row.bus_to)]] = row[direction]
    end
    @assert !any(isnan.(flows))
    return flows
end

@testset "SINGLE PERIOD power flows evaluation: ABA, PTDF, VirtualPTDF" begin
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # get reference values: flows and angles
    data = PowerFlowData(DCPowerFlow(), sys; correct_bustypes = true)
    power_injection =
        deepcopy(data.bus_activepower_injection - data.bus_activepower_withdrawals)
    matrix_data = deepcopy(data.power_network_matrix.K)       # LU factorization of ABA
    aux_network_matrix = deepcopy(data.aux_network_matrix)    # BA matrix

    valid_ix = setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
    ref_bus_angles = deepcopy(data.bus_angles)
    ref_flow_values = deepcopy(data.branch_activepower_flow_from_to)

    ref_bus_angles[valid_ix] = matrix_data \ power_injection[valid_ix]
    ref_flow_values = transpose(aux_network_matrix.data) * ref_bus_angles
    # ref_flow_values[i] gives the flow in the branch with arc data.aux_network_matrix.axes[2][i]
    # ref_bus_angles[i] gives the angle at bus with bus number data.power_network_matrix.axes[1][i]

    # CASE 1: ABA and BA matrices
    solved_data_ABA = solve_powerflow(DCPowerFlow(), sys; correct_bustypes = true)
    @test isapprox(
        flows_from_dataframe(solved_data_ABA, data.branch_lookup, :P_from_to),
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        flows_from_dataframe(solved_data_ABA, data.branch_lookup, :P_to_from),
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_ABA["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)

    # CASE 2: PTDF and ABA MATRICES
    solved_data_PTDF = solve_powerflow(PTDFDCPowerFlow(), sys; correct_bustypes = true)
    @test isapprox(
        flows_from_dataframe(solved_data_PTDF, data.branch_lookup, :P_from_to),
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        flows_from_dataframe(solved_data_PTDF, data.branch_lookup, :P_to_from),
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_PTDF["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)

    # CASE 3: VirtualPTDF and ABA MATRICES
    solved_data_vPTDF = solve_powerflow(vPTDFDCPowerFlow(), sys; correct_bustypes = true)
    @test isapprox(
        flows_from_dataframe(solved_data_vPTDF, data.branch_lookup, :P_from_to),
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        flows_from_dataframe(solved_data_vPTDF, data.branch_lookup, :P_to_from),
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_vPTDF["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)
end
