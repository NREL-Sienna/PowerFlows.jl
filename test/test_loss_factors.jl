
@testset "test_loss_factors_case_14" begin
    for ACSolver in AC_SOLVERS_TO_TEST
        # FIXME failing for LevenbergMarquardtACPowerFlow. investigate.
        @testset "AC Solver: $(ACSolver)" begin
            sys = build_system(PSITestSystems, "c_sys14"; add_forecasts = false)

            pf = ACPowerFlow(ACSolver)
            pf_lf = ACPowerFlow(ACSolver; calculate_loss_factors = true)
            time_steps = 24
            data_loss_factors =
                PowerFlowData(pf_lf, sys; time_steps = time_steps)
            data_brute_force =
                PowerFlowData(pf, sys; time_steps = time_steps)

            # allocate timeseries data from csv
            prepare_ts_data!(data_loss_factors, time_steps)
            prepare_ts_data!(data_brute_force, time_steps)

            # get power flows with NR KLU method and write results
            solve_powerflow!(data_loss_factors; pf = pf_lf)

            # get loss factors using brute force approach (sequential power flow evaluations for each bus)
            bf_loss_factors = penalty_factors_brute_force(data_brute_force, pf)

            # confirm that loss factors match for the Jacobian-based and brute force approaches
            @test all(
                isapprox.(
                    bf_loss_factors,
                    data_loss_factors.loss_factors,
                    atol = 1e-4,
                    rtol = 0,
                ),
            )

            # get power flow results without loss factors
            solve_powerflow!(data_brute_force; pf = pf)
            @test isnothing(data_brute_force.loss_factors)
        end
    end
end
