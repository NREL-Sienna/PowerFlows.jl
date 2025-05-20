@testset "NewtonRaphsonACPowerFlow" begin
    # test NR kwargs.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    nr_pf = ACPowerFlow{NewtonRaphsonACPowerFlow}()
    @test_logs (:info, r".*NewtonRaphsonACPowerFlow solver converged"
    ) match_mode = :any PF.solve_powerflow(nr_pf, sys; maxIterations = 50,
        tol = 1e-10, refinement_threshold = 0.01, refinement_eps = 1e-7)
end

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

# TODO: can I create an input such that newton doesn't converge, but iterative_refinement does?
# need jacobian to be ill-conditioned...

function bad_x0!(sys::PSY.System)
    for comp in get_components(PSY.PowerLoad, sys)
        set_angle!(PSY.get_bus(comp), 1.0)
    end
end

@testset "dc fallback" begin
    pf = ACPowerFlow{NewtonRaphsonACPowerFlow}(; robust_power_flow = true)
    # test that _dc_powerflow_fallback! solves correctly.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    sys2 = deepcopy(sys)
    data = PowerFlowData(pf, sys2)
    PF._dc_powerflow_fallback!(data, 1)
    ABA_angles = data.bus_angles[data.valid_ix, 1]
    p_inj =
        data.bus_activepower_injection[data.valid_ix, 1] -
        data.bus_activepower_withdrawals[data.valid_ix, 1]
    @test data.aux_network_matrix.data * ABA_angles ≈ p_inj

    # check behavior of improved_x0 via creating bogus awful starting point.
    sys3 = deepcopy(sys)
    bad_x0!(sys3)
    data3 = PowerFlowData(pf, sys3)
    x0 = PF.calculate_x0(data3, 1)
    residual = PF.ACPowerFlowResidual(data3, 1)
    residual(x0, 1)
    residualSize = norm(residual.Rv, 1)
    newx0 = deepcopy(x0)
    PF.improve_x0!(newx0, data, 1, residual)
    residual(newx0, 1)
    newResidualSize = norm(residual.Rv, 1)
    @test x0 !== newx0
    @test newResidualSize < residualSize
    # TODO: case with bad residual where DC powerflow doesn't yield improvement?

    # check that it does the DC fallback.
    # _initialize_bus_data! corrects the voltages to be "reasonable," between 0.8 and 1.2
    sys4 = deepcopy(sys)
    bad_x0!(sys4)
    improvement_regex = r".*DC powerflow fallback yields better x0"
    @test_logs (:info, improvement_regex) match_mode = :any PF.solve_powerflow(pf, sys4)
    pf_no_dc = ACPowerFlow{NewtonRaphsonACPowerFlow}(; robust_power_flow = false)
    sys5 = deepcopy(sys)
    bad_x0!(sys5)
    @test_logs (:debug, "skipping DC powerflow fallback") match_mode = :any min_level =
        Logging.Debug PF.solve_powerflow(pf_no_dc, sys5)
end

@testset "large residual warning" begin
    # a system where bus numbers aren't 1, 2,...is there a smaller one with this property?
    sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys")
    for i in [1, 35, 52, 57, 43, 66, 49, 68, 71, 25, 69, 58, 3, 73]
        data = PowerFlowData(ACPowerFlow(), sys)
        # First, write solution to data. Then set magnitude of a random-ish bus to a huge number
        # and try to solve again: the "large residual warning" should be about that bus.
        solve_powerflow!(data)
        bus = collect(PSY.get_components(PSY.ACBus, sys))[i]
        bus_no = bus.number
        bus_ix = PowerFlows.get_bus_lookup(data)[bus_no]
        data.bus_magnitude[bus_ix] = 100.0
        @test_logs (:warn, Regex(".*Largest residual at bus $(bus_no).*")
        ) match_mode = :any solve_powerflow!(data)
    end
end

# It happened: our more efficient "update J" means this test no longer finds 
# a point where J is singular.
#=
@testset "singular Jacobian trust region" begin
    # NewtonRaphsonACPowerFlow fails to converge on this system.
    # Empirically found: this system happens to have singular J's right next to the solution,
    # and if we call solve_powerflow! repeatedly, TrustRegion happens to pick such a point.
    # (May break if trust region is changed. TODO eventually: find a better test case.)
    pf = ACPowerFlow{TrustRegionACPowerFlow}()
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys")
    data = PowerFlowData(pf, sys)
    @test_logs (:warn, Regex(".*Jacobian is singular.*")
    ) match_mode = :any for _ in 1:20
        solve_powerflow!(data; pf = pf)
    end
end=#
