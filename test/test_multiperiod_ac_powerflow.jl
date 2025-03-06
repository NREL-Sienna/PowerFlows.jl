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
@testset "MULTI-PERIOD power flows evaluation: NR" for ACSolver in AC_SOLVERS_TO_TEST
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # create structure for multi-period case
    pf = ACPowerFlow(ACSolver)
    time_steps = 24
    data = PowerFlowData(pf, sys; time_steps = time_steps)

    # allocate timeseries data from csv
    prepare_ts_data!(data, time_steps)

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

@testset "MULTI-PERIOD power flows evaluation: compare results for different solvers" begin
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # create structure for multi-period case
    pf_lu = ACPowerFlow(PowerFlows.LUACPowerFlow)
    pf_test = ACPowerFlow(NewtonRaphsonACPowerFlow)

    time_steps = 24
    data_lu = PowerFlowData(pf_lu, sys; time_steps = time_steps)
    data_test = PowerFlowData(pf_test, sys; time_steps = time_steps)

    # allocate timeseries data from csv
    prepare_ts_data!(data_lu, time_steps)
    prepare_ts_data!(data_test, time_steps)

    # get power flows with NR KLU method and write results
    solve_powerflow!(data_lu; pf = pf_lu)
    solve_powerflow!(data_test; pf = pf_test)

    # check results
    @test isapprox(data_lu.bus_magnitude, data_test.bus_magnitude, atol = 1e-9)
    @test isapprox(data_lu.bus_angles, data_test.bus_angles, atol = 1e-9)
    @test isapprox(
        data_lu.branch_activepower_flow_from_to,
        data_test.branch_activepower_flow_from_to,
        atol = 1e-9,
    )
    @test isapprox(
        data_lu.branch_activepower_flow_to_from,
        data_test.branch_activepower_flow_to_from,
        atol = 1e-9,
    )
    @test isapprox(
        data_lu.branch_reactivepower_flow_from_to,
        data_test.branch_reactivepower_flow_from_to,
        atol = 1e-9,
    )
    @test isapprox(
        data_lu.branch_reactivepower_flow_to_from,
        data_test.branch_reactivepower_flow_to_from,
        atol = 1e-9,
    )
end

@testset "test_loss_factors_case_14" for ACSolver in AC_SOLVERS_TO_TEST
    sys = build_system(PSITestSystems, "c_sys14"; add_forecasts = false)

    pf = ACPowerFlow(ACSolver)
    time_steps = 24
    data_loss_factors =
        PowerFlowData(pf, sys; time_steps = time_steps, calculate_loss_factors = true)
    data_brute_force =
        PowerFlowData(pf, sys; time_steps = time_steps, calculate_loss_factors = false)

    # allocate timeseries data from csv
    prepare_ts_data!(data_loss_factors, time_steps)
    prepare_ts_data!(data_brute_force, time_steps)

    # get power flows with NR KLU method and write results
    solve_powerflow!(data_loss_factors; pf = pf)

    # get loss factors using brute force approach (sequential power flow evaluations for each bus)
    bf_loss_factors = penalty_factors_brute_force(data_brute_force, pf)

    # confirm that loss factors match for the Jacobian-based and brute force approaches
    @test isapprox(bf_loss_factors, data_loss_factors.loss_factors, atol = 1e-4, rtol = 0)

    # get power flow results without loss factors
    solve_powerflow!(data_brute_force; pf = pf)
    @test isnothing(data_brute_force.loss_factors)
end
