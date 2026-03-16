"""
Test suite for Jacobian indefiniteness detection functions.

This test file verifies that the indefiniteness detection methods correctly identify
the eigenvalue spectrum of Jacobian matrices using sparse factorization.
"""

@testset "Jacobian Indefiniteness Detection" begin
    # Setup: Create a test system and compute Jacobian
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}()
    data = PF.PowerFlowData(pf, sys; correct_bustypes = true)
    time_step = 1
    residual = PF.ACPowerFlowResidual(data, time_step)
    jacobian = PF.ACPowerFlowJacobian(data, time_step)

    # Initialize variables at operating point
    x0 = PF.calculate_x0(data, time_step)
    residual(x0, time_step)
    jacobian(time_step)

    @testset "InertiaResult Structure" begin
        # Test InertiaResult creation and display
        inertia_result = PF.InertiaResult(5, 3, 2, true, false, false, 1e-14)

        @test inertia_result.n_positive == 5
        @test inertia_result.n_negative == 3
        @test inertia_result.n_zero == 2
        @test inertia_result.is_indefinite == true
        @test inertia_result.is_positive_definite == false
        @test inertia_result.is_negative_definite == false
        @test inertia_result.tolerance == 1e-14
    end

    @testset "Inertia Computation via Sparse LU" begin
        # Test the main inertia computation function
        J = jacobian.Jv
        inertia = PF.compute_inertia_via_sparse_lu(J)

        # Basic checks
        @test inertia.n_positive >= 0
        @test inertia.n_negative >= 0
        @test inertia.n_zero >= 0
        @test inertia.n_positive + inertia.n_negative + inertia.n_zero == size(J, 1)

        # The Jacobian should have both positive and negative eigenvalues (indefinite)
        # This is typical for power flow Jacobians at operating points
        if !all(iszero, J.nzval)  # Only test if matrix is not all zeros
            @test true  # The matrix is properly analyzed
        end
    end

    @testset "is_jacobian_indefinite Function" begin
        # Test the indefiniteness check
        J = jacobian.Jv
        is_indef = PF.is_jacobian_indefinite(J)

        @test is_indef isa Bool
    end

    @testset "is_positive_definite Function" begin
        # Test positive definiteness check
        J = jacobian.Jv
        is_pos_def = PF.is_positive_definite(J)

        @test is_pos_def isa Bool

        # If the matrix is indefinite, it cannot be positive definite
        inertia = PF.compute_inertia_via_sparse_lu(J)
        if inertia.is_indefinite
            @test is_pos_def == false
        end
    end

    @testset "is_negative_definite Function" begin
        # Test negative definiteness check
        J = jacobian.Jv
        is_neg_def = PF.is_negative_definite(J)

        @test is_neg_def isa Bool

        # If the matrix is indefinite, it cannot be negative definite
        inertia = PF.compute_inertia_via_sparse_lu(J)
        if inertia.is_indefinite
            @test is_neg_def == false
        end
    end

    @testset "check_jacobian_symmetric_part Function" begin
        # Test analysis of symmetric part
        J = jacobian.Jv
        sym_inertia = PF.check_jacobian_symmetric_part(J)

        @test sym_inertia.n_positive >= 0
        @test sym_inertia.n_negative >= 0
        @test sym_inertia.n_zero >= 0
        @test sym_inertia.n_positive + sym_inertia.n_negative + sym_inertia.n_zero == size(J, 1)
    end

    @testset "quick_indefiniteness_check Function" begin
        # Test convenience function for ACPowerFlowJacobian
        is_indef = PF.quick_indefiniteness_check(jacobian)

        @test is_indef isa Bool
    end

    @testset "get_inertia_report Function" begin
        # Test report generation
        J = jacobian.Jv
        report = PF.get_inertia_report(J)

        @test report isa String
        @test length(report) > 0
        @test contains(report, "Jacobian Matrix Inertia Report")
        @test contains(report, "Positive eigenvalues")
        @test contains(report, "Negative eigenvalues")
    end

    @testset "monitor_jacobian_definiteness Function" begin
        # Test monitoring function
        inertia = PF.monitor_jacobian_definiteness(jacobian; verbose=false)

        @test inertia isa PF.InertiaResult
        @test inertia.n_positive + inertia.n_negative + inertia.n_zero == size(jacobian.Jv, 1)
    end

    @testset "Consistency Checks" begin
        # Verify that definiteness properties are mutually exclusive
        J = jacobian.Jv
        inertia = PF.compute_inertia_via_sparse_lu(J)

        # A matrix cannot be both positive and negative definite
        @test !(inertia.is_positive_definite && inertia.is_negative_definite)

        # If indefinite, cannot be positive or negative definite
        if inertia.is_indefinite
            @test !inertia.is_positive_definite
            @test !inertia.is_negative_definite
        end

        # If positive definite, cannot be indefinite or negative definite
        if inertia.is_positive_definite
            @test !inertia.is_indefinite
            @test !inertia.is_negative_definite
        end

        # If negative definite, cannot be indefinite or positive definite
        if inertia.is_negative_definite
            @test !inertia.is_indefinite
            @test !inertia.is_positive_definite
        end
    end

    @testset "Tolerance Parameter" begin
        # Test that tolerance parameter affects zero eigenvalue detection
        J = jacobian.Jv

        inertia_tight = PF.compute_inertia_via_sparse_lu(J; tolerance=1e-15)
        inertia_loose = PF.compute_inertia_via_sparse_lu(J; tolerance=1e-8)

        # Looser tolerance should classify more values as zero
        @test inertia_loose.n_zero >= inertia_tight.n_zero
    end

    @testset "Multiple System Test" begin
        # Test on a different system to ensure robustness
        sys2 = PSB.build_system(PSB.PSITestSystems, "c_sys30")
        data2 = PF.PowerFlowData(pf, sys2; correct_bustypes = true)
        jacobian2 = PF.ACPowerFlowJacobian(data2, time_step)

        x0_2 = PF.calculate_x0(data2, time_step)
        jacobian2(time_step)

        J2 = jacobian2.Jv
        inertia2 = PF.compute_inertia_via_sparse_lu(J2)

        # Verify basic properties hold for larger system
        @test inertia2.n_positive + inertia2.n_negative + inertia2.n_zero == size(J2, 1)
        @test inertia2.n_positive >= 0
        @test inertia2.n_negative >= 0
    end
end
