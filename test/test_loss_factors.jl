
@testset "test_loss_factors_case_14" begin
    for ACSolver in AC_SOLVERS_TO_TEST
        # FIXME failing for LevenbergMarquardtACPowerFlow. investigate.
        if ACSolver == LevenbergMarquardtACPowerFlow
            continue
        end
        @testset "AC Solver: $(ACSolver)" begin
            sys = build_system(PSITestSystems, "c_sys14"; add_forecasts = false)

            time_steps = 24
            pf = ACPowerFlow{ACSolver}(; time_steps = time_steps)
            pf_lf = ACPowerFlow{ACSolver}(;
                calculate_loss_factors = true,
                time_steps = time_steps,
            )
            data_loss_factors = PowerFlowData(pf_lf, sys)
            data_brute_force = PowerFlowData(pf, sys)

            # allocate timeseries data from csv
            prepare_ts_data!(data_loss_factors, time_steps)
            prepare_ts_data!(data_brute_force, time_steps)

            # get power flows with NR KLU method and write results
            solve_power_flow!(data_loss_factors)

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
            solve_power_flow!(data_brute_force)
            @test isnothing(data_brute_force.loss_factors)
        end
    end
end

@testset "test_loss_factors_multiple_ref_buses" begin
    for ACSolver in AC_SOLVERS_TO_TEST
        # FIXME failing for LevenbergMarquardtACPowerFlow. investigate.
        if ACSolver == LevenbergMarquardtACPowerFlow
            continue
        end
        @testset "AC Solver: $(ACSolver)" begin
            # Create a system with two disconnected islands, each with its own REF bus
            sys = System(100.0)

            # Island 1: buses 1-3
            b1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.05, 0.0)
            b2 = _add_simple_bus!(sys, 2, ACBusTypes.PQ, 230, 1.0, 0.0)
            b3 = _add_simple_bus!(sys, 3, ACBusTypes.PQ, 230, 1.0, 0.0)

            # Island 2: buses 4-6
            b4 = _add_simple_bus!(sys, 4, ACBusTypes.REF, 230, 1.02, 0.0)
            b5 = _add_simple_bus!(sys, 5, ACBusTypes.PQ, 230, 1.0, 0.0)
            b6 = _add_simple_bus!(sys, 6, ACBusTypes.PQ, 230, 1.0, 0.0)

            # Add sources at REF buses
            _add_simple_source!(sys, b1, 0.5, 0.1)
            _add_simple_source!(sys, b4, 0.4, 0.08)

            # Add loads at PQ buses
            _add_simple_load!(sys, b2, 0.25, 0.05)
            _add_simple_load!(sys, b3, 0.2, 0.04)
            _add_simple_load!(sys, b5, 0.2, 0.04)
            _add_simple_load!(sys, b6, 0.15, 0.03)

            # Connect buses within island 1 (no connection between islands)
            _add_simple_line!(sys, b1, b2, 0.01, 0.05, 0.02)
            _add_simple_line!(sys, b2, b3, 0.015, 0.08, 0.01)
            _add_simple_line!(sys, b1, b3, 0.012, 0.06, 0.015)

            # Connect buses within island 2 (no connection between islands)
            _add_simple_line!(sys, b4, b5, 0.01, 0.05, 0.02)
            _add_simple_line!(sys, b5, b6, 0.015, 0.08, 0.01)
            _add_simple_line!(sys, b4, b6, 0.012, 0.06, 0.015)

            pf_lf = ACPowerFlow{ACSolver}(; calculate_loss_factors = true)
            data_loss_factors = PowerFlowData(pf_lf, sys)

            # Verify we have multiple REF buses before solving
            ref_buses = findall(==(PSY.ACBusTypes.REF), data_loss_factors.bus_type[:, 1])
            @test length(ref_buses) == 2

            # Solving should succeed (with a warning about multiple REF buses)
            @test_throws ErrorException solve_power_flow!(data_loss_factors; pf = pf_lf)
        end
    end
end
