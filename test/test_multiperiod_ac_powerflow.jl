function prepare_ts_data!(data::PowerFlowData, time_steps::Int64 = 24)
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
    # allocate data from csv
    injs = Matrix(injections)
    withs = Matrix(withdrawals)

    data.bus_activepower_injection .= deepcopy(injs[:, 1:time_steps])
    data.bus_activepower_withdrawals .= deepcopy(withs[:, 1:time_steps])
    return nothing
end

# work in progress
@testset "MULTI-PERIOD power flows evaluation: NR" for ACSolver in
                                                       (
    KLUACPowerFlow,
    PowerFlows.LUACPowerFlow,
    HybridACPowerFlow,
)
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # create structure for multi-period case
    pf = ACPowerFlow(ACSolver)
    time_steps = 24
    data = PowerFlowData(pf, sys; time_steps = time_steps)

    # allocate timeseries data from csv
    prepare_ts_data!(data, time_steps)

    # get power flows with NR KLU method and write results
    solve_powerflow!(data; pf = pf, enable_progress_bar = false)

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
    KLUACPowerFlow,
    PowerFlows.LUACPowerFlow,
    HybridACPowerFlow,
)
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # create structure for multi-period case
    pf_klu = ACPowerFlow(KLUACPowerFlow)
    pf_test = ACPowerFlow(ACSolver)

    time_steps = 24
    data_klu = PowerFlowData(pf_klu, sys; time_steps = time_steps)
    data_test = PowerFlowData(pf_test, sys; time_steps = time_steps)

    # allocate timeseries data from csv
    prepare_ts_data!(data_klu, time_steps)
    prepare_ts_data!(data_test, time_steps)

    # get power flows with NR KLU method and write results
    solve_powerflow!(data_klu; pf = pf_klu, enable_progress_bar = false)
    solve_powerflow!(data_test; pf = pf_test, enable_progress_bar = false)

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

@testset "test_loss_factors_case_14" for ACSolver in (KLUACPowerFlow, HybridACPowerFlow)
    sys = build_system(PSITestSystems, "c_sys14"; add_forecasts = false)

    pf_klu = ACPowerFlow(ACSolver; calc_loss_factors = true)
    pf_no_factors = ACPowerFlow(ACSolver)  # calc_loss_factors = false by default

    time_steps = 24
    data_klu = PowerFlowData(pf_klu, sys; time_steps = time_steps)
    data_no_factors = PowerFlowData(pf_no_factors, sys; time_steps = time_steps)

    # allocate timeseries data from csv
    prepare_ts_data!(data_klu, time_steps)
    prepare_ts_data!(data_no_factors, time_steps)

    # get power flows with NR KLU method and write results
    solve_powerflow!(data_klu; pf = pf_klu, enable_progress_bar = false)

    # get loss factors using brute force approach (sequential power flow evaluations for each bus)
    bf_loss_factors =
        PowerFlows.penalty_factors_brute_force(data_klu; enable_progress_bar = false)

    # confirm that loss factors match for the Jacobian-based and brute force approaches
    @test isapprox(bf_loss_factors, data_klu.loss_factors, atol = 1e-5, rtol = 0)

    # get power flow results without loss factors
    solve_powerflow!(data_no_factors; pf = pf_no_factors, enable_progress_bar = false)
    @test isnothing(data_no_factors.loss_factors)
end
