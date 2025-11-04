@testset "NewtonRaphsonACPowerFlow kwargs" begin
    # test NR kwargs.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    nr_pf = ACPowerFlow{NewtonRaphsonACPowerFlow}()
    @test_logs (:info, r".*NewtonRaphsonACPowerFlow solver converged"
    ) match_mode = :any PF.solve_powerflow(nr_pf, sys; maxIterations = 50,
        tol = 1e-10, refinement_threshold = 0.01, refinement_eps = 1e-7)
end

@testset "TrustRegionACPowerFlow kwargs" begin
    # test trust region kwargs.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    tr_pf = ACPowerFlow{TrustRegionACPowerFlow}()
    @test_logs (:info, r".*TrustRegionACPowerFlow solver converged"
    ) match_mode = :any PF.solve_powerflow(tr_pf, sys; eta = 1e-5,
        tol = 1e-10, factor = 1.1, maxIterations = 50)
end

function bad_x0!(sys::PSY.System)
    for comp in get_components(PSY.PowerLoad, sys)
        set_angle!(PSY.get_bus(comp), 1.0)
    end
end

@testset "TrustRegionACPowerFlow behavior" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    tr_pf = ACPowerFlow{TrustRegionACPowerFlow}()

    # Small trust region size => Cauchy or dogleg step
    @test_logs (:debug, r"(Dogleg step selected|Cauchy step selected)") match_mode = :any min_level =
        Logging.Debug PF.solve_powerflow(
        tr_pf,
        sys;
        factor = 0.01,
        maxIterations = 1,
    )

    # Large trust region size => Newton-Raphson step
    @test_logs (:debug, r"Newton-Raphson step selected.*") match_mode = :any min_level =
        Logging.Debug PF.solve_powerflow(
        tr_pf,
        sys;
        factor = 10.0,
        maxIterations = 1,
    )

    # Large eta => step rejected
    @test_logs (:debug, r"Step rejected.*") match_mode = :any min_level = Logging.Debug PF.solve_powerflow(
        tr_pf,
        sys;
        eta = 2.0,
        maxIterations = 1,
    )

    # Small eta => step accepted
    @test_logs (:debug, r"Step accepted.*") match_mode = :any min_level = Logging.Debug PF.solve_powerflow(
        tr_pf,
        sys;
        eta = 1e-6,
        maxIterations = 1,
    )
end

@testset "dc fallback" begin
    dc_pf = ACPowerFlow{NewtonRaphsonACPowerFlow}(;
        robust_power_flow = true,
        enhanced_flat_start = false,
    )
    no_dc_pf = ACPowerFlow{NewtonRaphsonACPowerFlow}(;
        robust_power_flow = false,
        enhanced_flat_start = false,
    )
    # test that _dc_powerflow_fallback! solves correctly.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    sys2 = deepcopy(sys)
    data = PowerFlowData(dc_pf, sys2)
    PF._dc_powerflow_fallback!(data, 1)
    valid_ix = PF.get_valid_ix(data)
    ABA_angles = data.bus_angles[valid_ix, 1]
    p_inj =
        data.bus_activepower_injection[valid_ix, 1] -
        data.bus_activepower_withdrawals[valid_ix, 1]
    @test data.aux_network_matrix.data * ABA_angles ≈ p_inj

    # check behavior of improved_x0 via creating bogus awful starting point.
    sys3 = deepcopy(sys)
    bad_x0!(sys3)
    data3 = PowerFlowData(no_dc_pf, sys3)
    x0 = PF.calculate_x0(data3, 1)
    residual = PF.ACPowerFlowResidual(data3, 1)
    residual(x0, 1)
    residualSize = norm(residual.Rv, 1)
    newx0 = deepcopy(x0)
    PF.dc_powerflow_start!(newx0, data, 1, residual)
    residual(newx0, 1)
    newResidualSize = norm(residual.Rv, 1)
    @test x0 !== newx0
    @test newResidualSize < residualSize
    # TODO: case with bad residual where DC powerflow doesn't yield improvement?

    # check that it does the DC fallback.
    # _initialize_bus_data! corrects the voltages to be "reasonable," between 0.8 and 1.2
    sys4 = deepcopy(sys)
    bad_x0!(sys4)
    improvement_regex = r".*DC powerflow fallback yields smaller residual.*"
    @test_logs (:info, improvement_regex) match_mode = :any PF.solve_powerflow(
        dc_pf,
        sys4,
    )
    sys5 = deepcopy(sys)
    bad_x0!(sys5)
    @test_logs (:debug, "skipping running DC powerflow fallback") match_mode = :any min_level =
        Logging.Debug PF.solve_powerflow(no_dc_pf, sys5)
end

@testset "large residual warning" begin
    # a system where bus numbers aren't 1, 2,...is there a smaller one with this property?
    sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys")
    for i in [1, 35, 52, 57, 43, 66, 49, 68, 71, 25, 69, 58, 3, 73]
        pf = ACPowerFlow(; enhanced_flat_start = false)
        data = PowerFlowData(pf, sys; correct_bustypes = true)
        # First, write solution to data. Then set magnitude of a random-ish bus to a huge number
        # and try to solve again: the "large residual warning" should be about that bus.
        solve_powerflow!(data)
        bus = collect(PSY.get_components(PSY.ACBus, sys))[i]
        bus_no = bus.number
        bus_ix = PowerFlows.get_bus_lookup(data)[bus_no]
        data.bus_magnitude[bus_ix] = 100.0
        @test_logs (:warn, Regex(".*Largest residual at bus $(bus_no).*")
        ) match_mode = :any solve_powerflow!(data; pf = pf)
    end
end
