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
    results = solve_powerflow!(data)

    # check results
    # for t in 1:length(data.timestep_map)
    #     res_t = solve_powerflow(pf, sys; time_step=t)  # does not work - ts data not set in sys
    #     flow_ft = res_t["flow_results"].P_from_to
    #     flow_tf = res_t["flow_results"].P_to_from
    #     ts_flow_ft = results[data.timestep_map[t]]["flow_results"].P_from_to
    #     ts_flow_tf = results[data.timestep_map[t]]["flow_results"].P_to_from
    #     @test isapprox(ts_flow_ft, flow_ft, atol = 1e-9)
    #     @test isapprox(ts_flow_tf, flow_tf, atol = 1e-9)
    # end
end
