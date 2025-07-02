@testset "NewtonRaphsonACPowerFlow" begin
    # test NR kwargs.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    nr_solver = NewtonRaphsonACPowerFlow(;
        max_iterations = 50,
        tolerance = 1e-10,
    )
    nr_pf = ACPowerFlow(nr_solver; check_reactive_power_limits = true)
    @test_logs (:info, r".*NewtonRaphsonACPowerFlow solver converged"
    ) match_mode = :any PF.solve_power_flow(sys, nr_pf)
end

@testset "TrustRegionACPowerFlow" begin
    # test trust region kwargs.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    tr_solver = TrustRegionACPowerFlow(;
        max_iterations = 50,
        tolerance = 1e-10,
        eta = 1e-5,
        trust_region_factor = 1.1,
    )
    tr_pf = ACPowerFlow(tr_solver)
    @test_logs (:info, r".*TrustRegionACPowerFlow solver converged"
    ) match_mode = :any PF.solve_power_flow(sys, tr_pf)
    # TODO better tests? i.e. more granularly compare behavior to expected, not just check end result.
    # could check behavior of delta, ie that delta is increased/decreased properly.
end

@testset "large residual warning" begin
    # a system where bus numbers aren't 1, 2,...is there a smaller one with this property?
    sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys")
    for i in [1, 35, 52, 57, 43, 66, 49, 68, 71, 25, 69, 58, 3, 73]
        pf = ACPowerFlow()
        data = PowerFlowData(pf, sys; correct_bustypes = true)
        # First, write solution to data. Then set magnitude of a random-ish bus to a huge number
        # and try to solve again: the "large residual warning" should be about that bus.
        PF.solve_power_flow_data!(data, pf)
        bus = collect(PSY.get_components(PSY.ACBus, sys))[i]
        bus_no = bus.number
        bus_ix = PowerFlows.get_bus_lookup(data)[bus_no]
        data.bus_magnitude[bus_ix] = 100.0
        @test_logs (:warn, Regex(".*Largest residual at bus $(bus_no).*")
        ) match_mode = :any PF.solve_power_flow_data!(data, pf)
    end
end

#=
@testset "singular Jacobian trust region" begin
    # NewtonRaphsonACPowerFlow fails to converge on this system.
    # Empirically found: this system happens to have singular J's right next to the solution,
    # and if we call solve_power_flow_data! repeatedly, TrustRegion happens to pick such a point.
    # (May break if trust region is changed. TODO eventually: find a better test case.)
    pf = ACPowerFlow{TrustRegionACPowerFlow}()
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys")
    data = PowerFlowData(pf, sys, correct_bustypes = true)
    @test_logs (:warn, Regex(".*Jacobian is singular.*")
    ) match_mode = :any for _ in 1:20
        solve_power_flow_data!(data, pf)
    end
end
=#
