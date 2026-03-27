"""
Test suite for Jacobian pivot-sign diagnostics (indefiniteness detection heuristic).

These tests verify that the LU-based pivot-sign analysis correctly counts positive,
negative, and near-zero pivots. Note: for non-symmetric matrices, pivot signs are
a heuristic and do not correspond to eigenvalue signs.
"""

@testset "Jacobian Indefiniteness Detection" begin
    # Setup: Create a test system and compute Jacobian
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}(; correct_bustypes = true)
    data = PF.PowerFlowData(pf, sys)
    time_step = 1
    residual = PF.ACPowerFlowResidual(data, time_step)
    jacobian = PF.ACPowerFlowJacobian(data,
        residual.bus_slack_participation_factors,
        residual.subnetworks,
        time_step,
    )

    # Initialize variables at operating point
    x0 = PF.calculate_x0(data, time_step)
    residual(x0, time_step)
    jacobian(time_step)

    @testset "PivotSignResult Structure" begin
        r = PF.PivotSignResult(5, 3, 2, 1e-14, true)

        @test r.n_positive == 5
        @test r.n_negative == 3
        @test r.n_zero == 2
        @test r.tolerance == 1e-14
        @test r.success == true
        @test PF.has_mixed_sign_pivots(r)
        @test !PF.all_pivots_positive(r)
        @test !PF.all_pivots_negative(r)

        # InertiaResult alias still works
        @test PF.InertiaResult === PF.PivotSignResult
    end

    @testset "Pivot-Sign Computation via Sparse LU" begin
        J = jacobian.Jv
        result = PF.compute_inertia_via_sparse_lu(J)

        # Factorization should succeed on a well-conditioned Jacobian
        @test result.success

        # Pivot counts must sum to matrix dimension
        @test result.n_positive >= 0
        @test result.n_negative >= 0
        @test result.n_zero >= 0
        @test result.n_positive + result.n_negative + result.n_zero == size(J, 1)

        # A well-conditioned power flow Jacobian should have no near-zero pivots
        # and should have both positive and negative pivots (mixed-sign)
        @test result.n_zero == 0
        @test PF.has_mixed_sign_pivots(result)
        @test result.n_positive > 0
        @test result.n_negative > 0
        @test !PF.all_pivots_positive(result)
        @test !PF.all_pivots_negative(result)
    end

    @testset "Known Constructed Matrices" begin
        # Diagonal matrices: for diagonal M, LU gives U = M, so pivot signs == diagonal signs

        # All-positive diagonal
        M_pd = sparse(Int32[1, 2, 3], Int32[1, 2, 3], [2.0, 3.0, 4.0], 3, 3)
        r_pd = PF.compute_inertia_via_sparse_lu(M_pd)
        @test r_pd.success
        @test r_pd.n_positive == 3
        @test r_pd.n_negative == 0
        @test r_pd.n_zero == 0
        @test PF.all_pivots_positive(r_pd)
        @test !PF.all_pivots_negative(r_pd)
        @test !PF.has_mixed_sign_pivots(r_pd)

        # All-negative diagonal
        M_nd = sparse(Int32[1, 2, 3], Int32[1, 2, 3], [-2.0, -3.0, -4.0], 3, 3)
        r_nd = PF.compute_inertia_via_sparse_lu(M_nd)
        @test r_nd.success
        @test r_nd.n_positive == 0
        @test r_nd.n_negative == 3
        @test r_nd.n_zero == 0
        @test !PF.all_pivots_positive(r_nd)
        @test PF.all_pivots_negative(r_nd)
        @test !PF.has_mixed_sign_pivots(r_nd)

        # Mixed-sign diagonal
        M_id = sparse(Int32[1, 2, 3], Int32[1, 2, 3], [2.0, -3.0, 4.0], 3, 3)
        r_id = PF.compute_inertia_via_sparse_lu(M_id)
        @test r_id.success
        @test r_id.n_positive == 2
        @test r_id.n_negative == 1
        @test r_id.n_zero == 0
        @test !PF.all_pivots_positive(r_id)
        @test !PF.all_pivots_negative(r_id)
        @test PF.has_mixed_sign_pivots(r_id)

        # Singular diagonal (has a zero row/column)
        M_sg = sparse(Int32[1, 2], Int32[1, 2], [2.0, 4.0], 3, 3)
        r_sg = PF.compute_inertia_via_sparse_lu(M_sg)
        @test r_sg.n_positive + r_sg.n_negative + r_sg.n_zero == 3
        @test !PF.all_pivots_positive(r_sg)
        @test !PF.all_pivots_negative(r_sg)
        if r_sg.success
            @test r_sg.n_zero >= 1
        end
    end

    @testset "Factorization Failure Invariant" begin
        # When factorization fails, success == false and counts sum to n
        failed = PF.PivotSignResult(0, 0, 5, 1e-14, false)
        @test !failed.success
        @test failed.n_positive + failed.n_negative + failed.n_zero == 5
        # Computed functions must return false when success == false
        @test !PF.has_mixed_sign_pivots(failed)
        @test !PF.all_pivots_positive(failed)
        @test !PF.all_pivots_negative(failed)
    end

    @testset "is_jacobian_indefinite Function" begin
        J = jacobian.Jv
        is_mixed = PF.is_jacobian_indefinite(J)
        @test is_mixed isa Bool
        # c_sys14 Jacobian has mixed-sign pivots
        @test is_mixed == true
    end

    @testset "is_positive_definite Function" begin
        J = jacobian.Jv
        is_pos = PF.is_positive_definite(J)
        @test is_pos isa Bool
        # c_sys14 Jacobian has mixed-sign pivots, so cannot be all-positive
        @test is_pos == false
    end

    @testset "is_negative_definite Function" begin
        J = jacobian.Jv
        is_neg = PF.is_negative_definite(J)
        @test is_neg isa Bool
        # c_sys14 Jacobian has mixed-sign pivots, so cannot be all-negative
        @test is_neg == false
    end

    @testset "check_jacobian_symmetric_part Function" begin
        J = jacobian.Jv
        sym_result = PF.check_jacobian_symmetric_part(J)

        @test sym_result.success
        @test sym_result.n_positive >= 0
        @test sym_result.n_negative >= 0
        @test sym_result.n_zero >= 0
        @test sym_result.n_positive + sym_result.n_negative + sym_result.n_zero ==
              size(J, 1)
    end

    @testset "quick_indefiniteness_check Function" begin
        is_mixed = PF.quick_indefiniteness_check(jacobian)
        @test is_mixed isa Bool
    end

    @testset "get_inertia_report Function" begin
        J = jacobian.Jv
        report = PF.get_inertia_report(J)

        @test report isa String
        @test length(report) > 0
        @test contains(report, "Pivot-Sign Diagnostic Report")
        @test contains(report, "Positive pivots")
        @test contains(report, "Negative pivots")
    end

    @testset "monitor_jacobian_definiteness Function" begin
        result = PF.monitor_jacobian_definiteness(jacobian)

        @test result isa PF.PivotSignResult
        @test result.success
        @test result.n_positive + result.n_negative + result.n_zero ==
              size(jacobian.Jv, 1)
    end

    @testset "Consistency Checks" begin
        J = jacobian.Jv
        result = PF.compute_inertia_via_sparse_lu(J)

        # Cannot have all-positive and all-negative simultaneously
        @test !(PF.all_pivots_positive(result) && PF.all_pivots_negative(result))

        # Mixed-sign implies not all-positive and not all-negative
        if PF.has_mixed_sign_pivots(result)
            @test !PF.all_pivots_positive(result)
            @test !PF.all_pivots_negative(result)
        end

        # All-positive implies no mixed-sign and no all-negative
        if PF.all_pivots_positive(result)
            @test !PF.has_mixed_sign_pivots(result)
            @test !PF.all_pivots_negative(result)
        end

        # All-negative implies no mixed-sign and no all-positive
        if PF.all_pivots_negative(result)
            @test !PF.has_mixed_sign_pivots(result)
            @test !PF.all_pivots_positive(result)
        end
    end

    @testset "Tolerance Parameter" begin
        J = jacobian.Jv

        result_tight = PF.compute_inertia_via_sparse_lu(J; tolerance = 1e-15)
        result_loose = PF.compute_inertia_via_sparse_lu(J; tolerance = 1e-8)

        # Looser tolerance should classify more pivots as near-zero
        @test result_loose.n_zero >= result_tight.n_zero
    end

    @testset "Multiple System Test" begin
        sys2 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
        data2 = PF.PowerFlowData(pf, sys2)
        residual2 = PF.ACPowerFlowResidual(data2, time_step)
        jacobian2 = PF.ACPowerFlowJacobian(data2,
            residual2.bus_slack_participation_factors,
            residual2.subnetworks,
            time_step,
        )

        x0_2 = PF.calculate_x0(data2, time_step)
        residual2(x0_2, time_step)
        jacobian2(time_step)

        J2 = jacobian2.Jv
        result2 = PF.compute_inertia_via_sparse_lu(J2)

        @test result2.success
        @test result2.n_positive + result2.n_negative + result2.n_zero == size(J2, 1)
        @test result2.n_positive >= 0
        @test result2.n_negative >= 0
    end
end

@testset "monitor_jacobian solver integration" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    @testset "NewtonRaphsonACPowerFlow with monitor_jacobian" begin
        pf = ACPowerFlow{NewtonRaphsonACPowerFlow}(;
            correct_bustypes = true,
            solver_settings = Dict{Symbol, Any}(:monitor_jacobian => true),
        )
        result = solve_power_flow(pf, sys)
        @test result !== nothing
    end

    @testset "TrustRegionACPowerFlow with monitor_jacobian" begin
        pf = ACPowerFlow{TrustRegionACPowerFlow}(;
            correct_bustypes = true,
            solver_settings = Dict{Symbol, Any}(:monitor_jacobian => true),
        )
        result = solve_power_flow(pf, sys)
        @test result !== nothing
    end

    @testset "LevenbergMarquardtACPowerFlow with monitor_jacobian" begin
        pf = ACPowerFlow{LevenbergMarquardtACPowerFlow}(;
            correct_bustypes = true,
            solver_settings = Dict{Symbol, Any}(:monitor_jacobian => true),
        )
        result = solve_power_flow(pf, sys)
        @test result !== nothing
    end

    @testset "monitor_jacobian defaults to false" begin
        pf = ACPowerFlow{NewtonRaphsonACPowerFlow}(; correct_bustypes = true)
        result = solve_power_flow(pf, sys)
        @test result !== nothing
    end
end
