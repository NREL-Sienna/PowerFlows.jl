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

@testset "large residual warning" begin
    # a system where bus numbers aren't 1, 2,...is there a smaller one with this property?
    sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys")
    data = PowerFlowData(ACPowerFlow(), sys)
    # First, write solution to data. Then set magnitude of a random bus to a huge number
    # and try to solve again: the "large residual warning" should be about that bus.
    solve_powerflow!(data)
    bus = rand(collect(PSY.get_components(PSY.ACBus, sys)))
    bus_no = bus.number
    bus_ix = PowerFlows.get_bus_lookup(data)[bus_no]
    data.bus_magnitude[bus_ix] = 100.0
    @test_logs (:warn, Regex(".*Largest residual at bus $(bus_no).*")
    ) match_mode = :any solve_powerflow!(data)
end
