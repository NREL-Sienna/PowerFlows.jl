@testset "Jacobian Definiteness (Arpack)" begin
    @testset "Dense vs matrix-free operator on c_sys14" begin
        sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
        pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}(; correct_bustypes = true)
        data = PF.PowerFlowData(pf, sys)
        time_step = 1
        residual = PF.ACPowerFlowResidual(data, time_step)
        jacobian = PF.ACPowerFlowJacobian(
            data,
            residual.bus_slack_participation_factors,
            residual.subnetworks,
            time_step,
        )
        x0 = PF.calculate_x0(data, time_step)
        residual(x0, time_step)
        jacobian(time_step)

        cache = PF.ReducedJacobianCache(data, time_step)

        # Dense path (small system, so this is the fallback)
        result_dense = PF._check_definiteness_dense(cache, jacobian.Jv)

        # Force operator path: build the operator and apply dense eigen on it
        # to compare against the dense J_R
        op = PF.build_reduced_jacobian_operator(cache, jacobian.Jv)
        n = op.n
        # Materialize operator as dense matrix by applying to unit vectors
        J_R_from_op = zeros(n, n)
        e = zeros(n)
        for i in 1:n
            e[i] = 1.0
            LinearAlgebra.mul!(view(J_R_from_op, :, i), op, e)
            e[i] = 0.0
        end

        # Eigenvalues of materialized operator should match dense path
        S = LinearAlgebra.Symmetric((J_R_from_op + J_R_from_op') ./ 2)
        λs = LinearAlgebra.eigvals(S)
        λ_min_op = first(λs)
        λ_max_op = last(λs)
        _, idx = findmin(abs, λs)
        λ_nearest_zero_op = λs[idx]

        @test result_dense.smallest_eigenvalue ≈ λ_min_op atol = 1e-8
        @test result_dense.largest_eigenvalue ≈ λ_max_op atol = 1e-8
        @test result_dense.nearest_zero_eigenvalue ≈ λ_nearest_zero_op atol = 1e-8
        @test result_dense.is_positive_definite == (λ_min_op > 0)
        @test result_dense.is_negative_definite == (λ_max_op < 0)
    end

    @testset "Arpack vs dense on ACTIVSg2000" begin
        sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
        pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}(; correct_bustypes = true)
        data = PF.PowerFlowData(pf, sys)
        time_step = 1
        residual = PF.ACPowerFlowResidual(data, time_step)
        jacobian = PF.ACPowerFlowJacobian(
            data,
            residual.bus_slack_participation_factors,
            residual.subnetworks,
            time_step,
        )
        x0 = PF.calculate_x0(data, time_step)
        residual(x0, time_step)
        jacobian(time_step)

        cache = PF.ReducedJacobianCache(data, time_step)
        @test cache.n_pq > 10  # sanity: system should have many PQ buses

        # Arpack (matrix-free) path
        result_arpack = PF.check_definiteness(cache, jacobian.Jv)

        # Dense reference: materialize J_R and compute all eigenvalues
        op = PF.build_reduced_jacobian_operator(cache, jacobian.Jv)
        n = op.n
        J_R_dense = zeros(n, n)
        e = zeros(n)
        for i in 1:n
            e[i] = 1.0
            LinearAlgebra.mul!(view(J_R_dense, :, i), op, e)
            e[i] = 0.0
        end
        S = LinearAlgebra.Symmetric((J_R_dense + J_R_dense') ./ 2)
        λs = LinearAlgebra.eigvals(S)
        λ_min_dense = first(λs)
        λ_max_dense = last(λs)
        _, idx = findmin(abs, λs)
        λ_nearest_zero_dense = λs[idx]

        # Arpack eigenvalues should be close to dense eigenvalues
        @test result_arpack.smallest_eigenvalue ≈ λ_min_dense rtol = 1e-4
        @test result_arpack.largest_eigenvalue ≈ λ_max_dense rtol = 1e-4
        @test result_arpack.nearest_zero_eigenvalue ≈ λ_nearest_zero_dense rtol = 1e-4

        # Definiteness classification must agree
        @test result_arpack.is_positive_definite == (λ_min_dense > 0)
        @test result_arpack.is_negative_definite == (λ_max_dense < 0)
    end

    @testset "DefinitenessResult show method" begin
        r = PF.DefinitenessResult(true, false, 0.5, 10.0, 0.5)
        buf = IOBuffer()
        show(buf, r)
        output = String(take!(buf))
        @test contains(output, "POSITIVE DEFINITE")
        @test contains(output, "λ_min=0.5")
        @test contains(output, "λ_max=10.0")
        @test contains(output, "λ_nearest_zero=0.5")
    end

    @testset "get_definiteness_report (c_sys14)" begin
        sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
        pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}(; correct_bustypes = true)
        data = PF.PowerFlowData(pf, sys)
        time_step = 1
        residual = PF.ACPowerFlowResidual(data, time_step)
        jacobian = PF.ACPowerFlowJacobian(
            data,
            residual.bus_slack_participation_factors,
            residual.subnetworks,
            time_step,
        )
        x0 = PF.calculate_x0(data, time_step)
        residual(x0, time_step)
        jacobian(time_step)

        report = PF.get_definiteness_report(jacobian, time_step)
        @test report isa String
        @test contains(report, "Reduced Jacobian Definiteness Report")
        @test contains(report, "matrix-free")
        @test contains(report, "λ_min")
        @test contains(report, "λ_max")
        @test contains(report, "λ_nearest_zero")
    end
end
