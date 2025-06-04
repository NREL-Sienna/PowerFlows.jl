@testset "MULTI-PERIOD power flows evaluation: ABA, PTDF, VirtualPTDF" begin
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    injections = CSV.read(
        joinpath(MAIN_DIR, "test", "test_data", "c_sys14_injections.csv"),
        DataFrame;
        header = 0,
    )
    withdrawals = CSV.read(
        joinpath(MAIN_DIR, "test", "test_data", "c_sys14_withdrawals.csv"),
        DataFrame;
        header = 0,
    )
    flows = CSV.read(
        joinpath(MAIN_DIR, "test", "test_data", "c_sys14_flows.csv"),
        DataFrame;
        header = 0,
    )
    angles = CSV.read(
        joinpath(MAIN_DIR, "test", "test_data", "c_sys14_angles.csv"),
        DataFrame;
        header = 0,
    )

    ##############################################################################

    # create structure for multi-period case
    time_steps = 24
    data_1 = PowerFlowData(DCPowerFlow(), sys; time_steps = time_steps, fix_bustypes = true)
    data_2 =
        PowerFlowData(PTDFDCPowerFlow(), sys; time_steps = time_steps, fix_bustypes = true)
    data_3 =
        PowerFlowData(vPTDFDCPowerFlow(), sys; time_steps = time_steps, fix_bustypes = true)

    # allocate data from csv
    injs = Matrix(injections)
    withs = Matrix(withdrawals)

    data_1.bus_activepower_injection .= deepcopy(injs)
    data_1.bus_activepower_withdrawals .= deepcopy(withs)

    data_2.bus_activepower_injection .= deepcopy(injs)
    data_2.bus_activepower_withdrawals .= deepcopy(withs)

    data_3.bus_activepower_injection .= deepcopy(injs)
    data_3.bus_activepower_withdrawals .= deepcopy(withs)

    # case 1: get power flows with ABA method and write results
    results_1 = solve_powerflow(data_1, sys)

    # case 2: get power flows PTDF method and write results
    results_2 = solve_powerflow(data_2, sys)

    # case 3: get power flows Virtual PTDF method and write results
    results_3 = solve_powerflow(data_3, sys)

    # check results

    # case 1
    for i in 1:length(data_1.timestep_map)
        net_flow = results_1[data_1.timestep_map[i]]["flow_results"].P_from_to
        net_flow_tf = results_1[data_1.timestep_map[i]]["flow_results"].P_to_from
        @test isapprox(net_flow, flows[:, i], atol = 1e-5)
        @test isapprox(net_flow_tf, -flows[:, i], atol = 1e-5)
    end

    # case 2
    for i in 1:length(data_1.timestep_map)
        net_flow = results_1[data_2.timestep_map[i]]["flow_results"].P_from_to
        net_flow_tf = results_1[data_2.timestep_map[i]]["flow_results"].P_to_from
        @test isapprox(net_flow, flows[:, i], atol = 1e-5)
        @test isapprox(net_flow_tf, -flows[:, i], atol = 1e-5)
    end

    # case 3
    for i in 1:length(data_1.timestep_map)
        net_flow = results_1[data_3.timestep_map[i]]["flow_results"].P_from_to
        net_flow_tf = results_1[data_3.timestep_map[i]]["flow_results"].P_to_from
        @test isapprox(net_flow, flows[:, i], atol = 1e-5)
        @test isapprox(net_flow_tf, -flows[:, i], atol = 1e-5)
    end
end
