sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
function run_pf(loads::Vector{Float64}; kwargs...)
    sys2 = deepcopy(sys)
    @assert(length(get_components(PSY.PowerLoad, sys2)) == length(loads))
    for (i, comp) in enumerate(get_components(PSY.PowerLoad, sys2))
        set_max_active_power!(comp, loads[i])
        set_active_power!(comp, loads[i])
    end
    return PF.solve_powerflow(
        ACPowerFlow{TrustRegionACPowerFlow}(),
        sys2;
        kwargs...
    )
end

@testset "trust region" begin
    @test_logs (:info, r".*TrustRegionACPowerFlow solver converged") match_mode = :any run_pf(
        [
        35.407000101,
        23.5000028166,
        50.0,
    ])
    # test trust region kwargs.
    @test_logs (:info, r".*TrustRegionACPowerFlow solver converged") match_mode = :any run_pf(
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
