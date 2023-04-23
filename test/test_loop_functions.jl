@testset "power flows evaluation: ABA, PTDF, VirtualPTDF" begin
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # get reference values: flows and angles
    data = PowerFlowData(DCPowerFlow(), sys)
    power_injection =
        deepcopy(data.bus_activepower_injection - data.bus_activepower_withdrawals)
    matrix_data = deepcopy(data.power_network_matrix.K)       # LU factorization of ABA
    aux_network_matrix = deepcopy(data.aux_network_matrix)    # BA matrix

    valid_ix = setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
    ref_bus_angles = deepcopy(data.bus_angles)
    ref_flow_values = deepcopy(data.branch_flow_values)

    ref_bus_angles[valid_ix] = matrix_data \ power_injection[valid_ix]
    ref_flow_values = transpose(aux_network_matrix.data) * ref_bus_angles

    # CASE 1: ABA and BA matrices 
    solved_data_ABA = solve_powerflow!(DCPowerFlow(), sys)
    @test isapprox(solved_data_ABA.branch_flow_values, ref_flow_values, atol = 1e-6)
    @test isapprox(solved_data_ABA.bus_angles, ref_bus_angles, atol = 1e-6)

    # CASE 2: PTDF and ABA MATRICES
    solved_data_PTDF = solve_powerflow!(PTDFDCPowerFlow(), sys)
    @test isapprox(solved_data_PTDF.branch_flow_values, ref_flow_values, atol = 1e-6)
    @test isapprox(solved_data_PTDF.bus_angles, ref_bus_angles, atol = 1e-6)

    # CASE 3: VirtualPTDF and ABA MATRICES
    solved_data_vPTDF1 = PowerFlowData(vPTDFDCPowerFlow(), sys)
    @test isapprox(solved_data_vPTDF1.branch_flow_values, ref_flow_values, atol = 1e-6)
    @test isapprox(solved_data_vPTDF1.bus_angles, ref_bus_angles, atol = 1e-6)

    # TODO CASE 3: VirtualPTDF and ABA MATRICES (partial number of lines)

end
