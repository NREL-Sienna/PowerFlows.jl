# Work in progress
@testset "MULTI-PERIOD power flows evaluation: NR" for ACSolver in
                                                       (
    NLSolveACPowerFlow,
    KLUACPowerFlow,
    PowerFlows.LUACPowerFlow,
    HyrbidACPowerFlow
)
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
    pf = ACPowerFlow(ACSolver)
    time_steps = 24
    data = PowerFlowData(pf, sys; time_steps = time_steps)

    # allocate data from csv
    injs = Matrix(injections)
    withs = Matrix(withdrawals)

    data.bus_activepower_injection .= deepcopy(injs)
    data.bus_activepower_withdrawals .= deepcopy(withs)

    # get power flows with NR KLU method and write results
    solve_powerflow!(data; pf = pf)

    # check results
    # for t in 1:length(data.timestep_map)
    #     res_t = solve_powerflow(pf, sys, t)  # does not work - ts data not set in sys
    #     flow_ft = res_t["flow_results"].P_from_to
    #     flow_tf = res_t["flow_results"].P_to_from
    #     ts_flow_ft = results[data.timestep_map[t]]["flow_results"].P_from_to
    #     ts_flow_tf = results[data.timestep_map[t]]["flow_results"].P_to_from
    #     @test isapprox(ts_flow_ft, flow_ft, atol = 1e-9)
    #     @test isapprox(ts_flow_tf, flow_tf, atol = 1e-9)
    # end
end

@testset "MULTI-PERIOD power flows evaluation: compare results for different solvers" for ACSolver in
                                                                                          (
    NLSolveACPowerFlow,
    PowerFlows.LUACPowerFlow,
    HyrbidACPowerFlow
)
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
    pf_klu = ACPowerFlow(KLUACPowerFlow)
    pf_test = ACPowerFlow(ACSolver)

    time_steps = 24
    data_klu = PowerFlowData(pf_klu, sys; time_steps = time_steps)
    data_test = PowerFlowData(pf_test, sys; time_steps = time_steps)

    # allocate data from csv
    injs = Matrix(injections)
    withs = Matrix(withdrawals)

    data_klu.bus_activepower_injection .= deepcopy(injs[:, 1:time_steps])
    data_klu.bus_activepower_withdrawals .= deepcopy(withs[:, 1:time_steps])

    data_test.bus_activepower_injection .= deepcopy(injs[:, 1:time_steps])
    data_test.bus_activepower_withdrawals .= deepcopy(withs[:, 1:time_steps])

    # get power flows with NR KLU method and write results
    solve_powerflow!(data_klu; pf = pf_klu)
    solve_powerflow!(data_test; pf = pf_test)

    # check results
    @test isapprox(data_klu.bus_magnitude, data_test.bus_magnitude, atol = 1e-9)
    @test isapprox(data_klu.bus_angles, data_test.bus_angles, atol = 1e-9)
    @test isapprox(
        data_klu.branch_activepower_flow_from_to,
        data_test.branch_activepower_flow_from_to,
        atol = 1e-9,
    )
    @test isapprox(
        data_klu.branch_activepower_flow_to_from,
        data_test.branch_activepower_flow_to_from,
        atol = 1e-9,
    )
    @test isapprox(
        data_klu.branch_reactivepower_flow_from_to,
        data_test.branch_reactivepower_flow_from_to,
        atol = 1e-9,
    )
    @test isapprox(
        data_klu.branch_reactivepower_flow_to_from,
        data_test.branch_reactivepower_flow_to_from,
        atol = 1e-9,
    )
end
