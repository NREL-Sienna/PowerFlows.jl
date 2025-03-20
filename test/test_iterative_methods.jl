@testset "TrustRegionACPowerFlow" begin
    # test trust region kwargs.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    tr_pf = ACPowerFlow{TrustRegionACPowerFlow}()
    @test_logs (:info, r".*TrustRegionACPowerFlow solver converged"
    ) match_mode = :any PF.solve_powerflow(tr_pf, sys; eta = 1e-5,
        tol = 1e-10, factor = 1.1, maxIterations = 50)
    # TODO better tests? i.e. more granularly compare behavior to expected, not just check end result.
    # could check behavior of delta, ie that delta is increased/decreased properly.
end

@testset "NewtonRaphsonACPowerFlow" begin
    # test NR kwargs.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    nr_pf = ACPowerFlow{NewtonRaphsonACPowerFlow}()
    @test_logs (:info, r".*NewtonRaphsonACPowerFlow solver converged"
    ) match_mode = :any PF.solve_powerflow(nr_pf, sys; maxIterations = 50,
        tol = 1e-10, refinement_threshold = 0.01, refinement_eps = 1e-7)
end
