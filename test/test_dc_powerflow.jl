@testset "SINGLE PERIOD power flows evaluation: ABA, PTDF, VirtualPTDF" begin
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # get indices
    buses = collect(PSY.get_components(PSY.ACBus, sys))
    idx = sortperm(buses; by = x -> PSY.get_number(x))

    # get sorted indices for branches
    branches = collect(PSY.get_components(PSY.get_available, PSY.ACBranch, sys))
    from_bus = [PSY.get_number(PSY.get_arc(x).from) for x in branches]
    idx2 = sortperm(from_bus)

    # get reference values: flows and angles
    data = PowerFlowData(DCPowerFlow(), sys; correct_bustypes = true)
    power_injection =
        deepcopy(data.bus_activepower_injection - data.bus_activepower_withdrawals)
    matrix_data = deepcopy(data.power_network_matrix.K)       # LU factorization of ABA
    aux_network_matrix = deepcopy(data.aux_network_matrix)    # BA matrix

    valid_ix = setdiff(
        1:length(power_injection),
        PNM.get_ref_bus_position(data.aux_network_matrix),
    )
    ref_bus_angles = deepcopy(data.bus_angles)
    ref_flow_values = deepcopy(data.branch_activepower_flow_from_to)

    ref_bus_angles[valid_ix] = matrix_data \ power_injection[valid_ix]
    ref_flow_values = transpose(aux_network_matrix.data) * ref_bus_angles

    # CASE 1: ABA and BA matrices
    solved_data_ABA = solve_powerflow(DCPowerFlow(), sys; correct_bustypes = true)
    @test isapprox(
        solved_data_ABA["1"]["flow_results"].P_from_to,
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        solved_data_ABA["1"]["flow_results"].P_to_from,
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_ABA["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)

    # CASE 2: PTDF and ABA MATRICES
    solved_data_PTDF = solve_powerflow(PTDFDCPowerFlow(), sys; correct_bustypes = true)
    @test isapprox(
        solved_data_PTDF["1"]["flow_results"].P_from_to,
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        solved_data_PTDF["1"]["flow_results"].P_to_from,
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_PTDF["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)

    # CASE 3: VirtualPTDF and ABA MATRICES
    solved_data_vPTDF = solve_powerflow(vPTDFDCPowerFlow(), sys; correct_bustypes = true)
    @test isapprox(
        solved_data_vPTDF["1"]["flow_results"].P_from_to,
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        solved_data_vPTDF["1"]["flow_results"].P_to_from,
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_vPTDF["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)
end
