@testset "Reduced Jacobian Definiteness Detection" begin
    @testset "DefinitenessResult show" begin
        r = PF.DefinitenessResult(true, false, 0.5, 10.0, 0.5)
        buf = IOBuffer()
        show(buf, r)
        output = String(take!(buf))
        @test contains(output, "POSITIVE DEFINITE")
        @test contains(output, "λ_min=0.5")
        @test contains(output, "λ_nearest_zero=0.5")
    end

    @testset "Reduced Jacobian from power flow (c_sys14)" begin
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

        # Check definiteness via the Jacobian wrapper
        result = PF.check_definiteness(jacobian, time_step)
        @test result isa PF.DefinitenessResult
        # A well-conditioned, converging power flow should have a PD reduced Jacobian
        @test result.is_positive_definite
        @test result.smallest_eigenvalue > 0
        @test result.largest_eigenvalue > 0
    end

    @testset "c_sys5 reduced Jacobian" begin
        sys2 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
        pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}(; correct_bustypes = true)
        data2 = PF.PowerFlowData(pf, sys2)
        time_step = 1
        residual2 = PF.ACPowerFlowResidual(data2, time_step)
        jacobian2 = PF.ACPowerFlowJacobian(data2,
            residual2.bus_slack_participation_factors,
            residual2.subnetworks,
            time_step,
        )
        x0_2 = PF.calculate_x0(data2, time_step)
        residual2(x0_2, time_step)
        jacobian2(time_step)

        result2 = PF.check_definiteness(jacobian2, time_step)
        @test result2 isa PF.DefinitenessResult
    end

    @testset "get_definiteness_report" begin
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

        report = PF.get_definiteness_report(jacobian, time_step)
        @test contains(report, "Reduced Jacobian")
        @test contains(report, "Definiteness")
        @test contains(report, "λ_min")
        @test contains(report, "λ_max")
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
