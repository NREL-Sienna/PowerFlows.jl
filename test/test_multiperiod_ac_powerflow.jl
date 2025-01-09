# Work in progress
@testset "MULTI-PERIOD power flows evaluation: NR" begin
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
    pf = ACPowerFlow()
    time_steps = 24
    data = PowerFlowData(pf, sys; time_steps = time_steps)

    # allocate data from csv
    injs = Matrix(injections)
    withs = Matrix(withdrawals)

    data.bus_activepower_injection .= deepcopy(injs)
    data.bus_activepower_withdrawals .= deepcopy(withs)

    # get power flows with NR KLU method and write results
    ts_converged, ts_x = solve_powerflow(pf, data, sys)

    # todo: implement write_results for multiperiod power flow
    # todo: test the result values

    # check results
    # for i in 1:length(data.timestep_map)
    #     net_flow = results[data.timestep_map[i]]["flow_results"].P_from_to
    #     net_flow_tf = results[data.timestep_map[i]]["flow_results"].P_to_from
    #     @test isapprox(net_flow, flows[:, i], atol = 1e-5)
    #     @test isapprox(net_flow_tf, -flows[:, i], atol = 1e-5)
    # end
end
