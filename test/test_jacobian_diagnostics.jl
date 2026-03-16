@testset "Jacobian Diagnostics" begin
    @testset "Well-conditioned matrix" begin
        Random.seed!(42)
        n = 50
        A = SparseMatrixCSC{Float64, Int32}(
            sparse(1.0 * LinearAlgebra.I, n, n) + 0.1 * sprandn(Float64, n, n, 0.1),
        )
        jd = PF.compute_jacobian_diagnostics(A)

        @test jd.n == n
        @test jd.numerical_rank == n
        @test jd.is_singular == false
        @test jd.is_ill_conditioned == false
        @test jd.rcond > 1e-3
        @test jd.condest < 1e6
        @test jd.min_pivot > 0.0
        @test jd.max_pivot > 0.0
    end

    @testset "Singular matrix" begin
        n = 10
        # Row 10 = Row 1, making the matrix rank-deficient
        A = SparseMatrixCSC{Float64, Int32}(sparse(1.0 * LinearAlgebra.I, n, n))
        A[n, :] = A[1, :]
        jd = PF.compute_jacobian_diagnostics(A)

        @test jd.is_singular == true
        @test jd.numerical_rank < n
        @test jd.rcond < 1e-12
    end

    @testset "Ill-conditioned matrix" begin
        # Build a non-diagonal ill-conditioned matrix so KLU pivoting can't trivialize it.
        # Use a tridiagonal with entries scaled to create large condition number.
        n = 20
        I_idx = Int32[]
        J_idx = Int32[]
        V = Float64[]
        for i in 1:n
            push!(I_idx, Int32(i)); push!(J_idx, Int32(i)); push!(V, 10.0^(-i + 1))
            if i < n
                push!(I_idx, Int32(i)); push!(J_idx, Int32(i + 1)); push!(V, 1e-10)
                push!(I_idx, Int32(i + 1)); push!(J_idx, Int32(i)); push!(V, 1e-10)
            end
        end
        A = SparseMatrixCSC{Float64, Int32}(sparse(I_idx, J_idx, V, n, n))
        jd = PF.compute_jacobian_diagnostics(A)

        @test jd.numerical_rank == n
        @test jd.is_singular == false
        # rcond should be very small for this ill-conditioned matrix
        @test jd.rcond < 1e-4
    end

    @testset "Cache-based vs standalone consistency" begin
        Random.seed!(123)
        n = 30
        A = SparseMatrixCSC{Float64, Int32}(
            sparse(1.0 * LinearAlgebra.I, n, n) + 0.05 * sprandn(Float64, n, n, 0.2),
        )

        # Standalone
        jd_standalone = PF.compute_jacobian_diagnostics(A)

        # Cache-based
        cache = PF.KLULinSolveCache(A)
        PF.full_factor!(cache, A)
        jd_cache = PF.compute_jacobian_diagnostics(cache, A)

        @test jd_standalone.rcond ≈ jd_cache.rcond
        @test jd_standalone.condest ≈ jd_cache.condest
        @test jd_standalone.rgrowth ≈ jd_cache.rgrowth
        @test jd_standalone.numerical_rank == jd_cache.numerical_rank
        @test jd_standalone.is_singular == jd_cache.is_singular
        @test jd_standalone.is_ill_conditioned == jd_cache.is_ill_conditioned
    end

    @testset "quick_jacobian_check" begin
        # Healthy matrix
        n = 20
        A_good = SparseMatrixCSC{Float64, Int32}(sparse(1.0 * LinearAlgebra.I, n, n))
        cache_good = PF.KLULinSolveCache(A_good)
        PF.full_factor!(cache_good, A_good)
        @test PF.quick_jacobian_check(cache_good) == true

        # Near-singular: use standalone check on the singular matrix
        # (compute_jacobian_diagnostics catches the SingularException)
        A_sing = SparseMatrixCSC{Float64, Int32}(sparse(1.0 * LinearAlgebra.I, n, n))
        A_sing[n, :] = A_sing[1, :]  # rank deficient
        jd_bad = PF.compute_jacobian_diagnostics(A_sing)
        @test jd_bad.is_singular == true
        @test jd_bad.is_ill_conditioned == true
    end

    @testset "Tolerance parameter" begin
        n = 10
        # Moderately ill-conditioned: rcond ~ 1e-8
        vals = [10.0^(-i * 0.9) for i in 0:(n - 1)]
        A = SparseMatrixCSC{Float64, Int32}(sparse(1:n, 1:n, vals, n, n))

        jd_tight = PF.compute_jacobian_diagnostics(A; tolerance = 1e-15)
        jd_loose = PF.compute_jacobian_diagnostics(A; tolerance = 1e-2)

        # Same rcond value regardless of tolerance
        @test jd_tight.rcond ≈ jd_loose.rcond
        # Loose tolerance should flag more matrices as ill-conditioned
        if jd_tight.is_ill_conditioned
            @test jd_loose.is_ill_conditioned == true
        end
    end

    @testset "get_diagnostics_report" begin
        n = 10
        A = SparseMatrixCSC{Float64, Int32}(sparse(1.0 * LinearAlgebra.I, n, n))
        jd = PF.compute_jacobian_diagnostics(A)
        report = PF.get_diagnostics_report(jd)

        @test report isa String
        @test contains(report, "Jacobian Diagnostics Report")
        @test contains(report, "Reciprocal condition number")
        @test contains(report, "Pivot growth factor")
        @test contains(report, "OK")
    end

    @testset "JacobianDiagnostics show method" begin
        jd = PF.JacobianDiagnostics(0.5, 2.0, 1.0, 10, 10, false, false, 0.1, 5.0)
        buf = IOBuffer()
        show(buf, jd)
        output = String(take!(buf))
        @test contains(output, "JacobianDiagnostics")
        @test contains(output, "OK")
    end

    @testset "Power flow Jacobian (c_sys14)" begin
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
        x0 = PF.calculate_x0(data, time_step)
        residual(x0, time_step)
        jacobian(time_step)

        J = jacobian.Jv
        jd = PF.compute_jacobian_diagnostics(J)

        @test jd.n == size(J, 1)
        @test jd.numerical_rank == size(J, 1)
        @test jd.is_singular == false
        @test jd.rcond > 0.0
        @test jd.condest > 1.0
        @test jd.min_pivot > 0.0

        # quick_jacobian_check on ACPowerFlowJacobian
        @test PF.quick_jacobian_check(jacobian) isa Bool
    end
end
