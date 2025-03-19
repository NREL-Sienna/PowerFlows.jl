sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
function run_pf(loads::Vector{Float64}; kwargs...)
    sys2 = deepcopy(sys)
    @assert(length(get_components(PSY.PowerLoad, sys2)) == length(loads))
    for (i, comp) in enumerate(get_components(PSY.PowerLoad, sys2))
        set_max_active_power!(comp, loads[i])
        set_active_power!(comp, loads[i])
    end
    return PF.solve_powerflow(
        ACPowerFlow{PF.NewtonRaphsonACPowerFlow}(),
        sys2;
        kwargs,
    )
end

@testset "trust region" begin
    # found by hand: all methods fail to converge => decrease loads, newton converges => increase loads.
    # repeat, until find a point where trust region method converges and newton doesn't.
    # (if these values stop working, could automate the above process via taking midpoint or centroid.)
    @test_logs (:info, r"converged.*TrustRegionNRMethod") match_mode = :any run_pf(
        [
        35.407000101,
        23.5000028166,
        50.0,
    ])
    # test trust region kwargs.
    @test_logs (:info, r"converged.*TrustRegionNRMethod") match_mode = :any run_pf(
        [
            35.407000101,
            23.5000028166,
            50.0,
        ]; maxIterations = 30, eta = 1e-4, factor = 1.0)
    # TODO better tests? i.e. more granularly compare behavior to expected, not just check end result.
    # could check behavior of delta, ie that delta is increased/decreased properly.
end

# TODO: can I create an input such that newton doesn't converge, but iterative_refinement does?
# need jacobian to be ill-conditioned...

@testset "dc fallback" begin
    pf = ACPowerFlow{NewtonRaphsonACPowerFlow}()
    # test that _dc_powerflow_fallback! solves correctly.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    sys2 = deepcopy(sys)
    data = PowerFlowData(pf, sys2)
    PF._dc_powerflow_fallback!(data, 1)
    ABA_angles = data.bus_angles[data.valid_ix, 1]
    p_inj = data.bus_activepower_injection[data.valid_ix, 1] 
            - data.bus_activepower_withdrawals[data.valid_ix, 1]
    @test data.aux_network_matrix.data * ABA_angles ≈ p_inj

    # check behavior of improved_x0 via creating bogus awful starting point.
    for (i, comp) in enumerate(get_components(PSY.PowerLoad, sys2))
        set_max_active_power!(comp, 100.0)
        set_active_power!(comp, 100.0)
    end
    data = PowerFlowData(pf, sys2)
    x0 = PF.calculate_x0(data, 1)
    residual = PF.ACPowerFlowResidual(data, 1)
    residual(x0, 1)
    residualSize = norm(residual.Rv, 1)
    newx0 = PF.improved_x0(data, 1, residual)
    newResidualSize = norm(residual.Rv, 1)
    @test x0 !== newx0
    @test newResidualSize < residualSize
    # TODO: case with bad residual where DC powerflow doesn't yield improvement?

    # check that it does the DC fallback.
    for (i, comp) in enumerate(get_components(PSY.PowerLoad, sys2))
        set_max_active_power!(comp, 100.0)
        set_active_power!(comp, 100.0)
    end
    improvement_regex = r".*DC powerflow fallback yields better x0"
    @test_logs (:info, improvement_regex) match_mode = :any PF.solve_powerflow(pf, sys2)
end