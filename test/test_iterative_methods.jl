sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
function run_pf(loads::Vector{Float64}; maxIterations::Int = 30)
    sys2 = deepcopy(sys)
    @assert(length(get_components(PSY.PowerLoad, sys2)) == length(loads))
    for (i, comp) in enumerate(get_components(PSY.PowerLoad, sys2))
        set_max_active_power!(comp, loads[i])
        set_active_power!(comp, loads[i])
    end
    return PF.solve_powerflow(
        ACPowerFlow{PF.HybridACPowerFlow}(),
        sys2;
        maxIterations = maxIterations,
    )
end

@testset "trust region" begin
    # found by hand: all methods fail to converge => decrease loads, newton converges => increase loads.
    # repeat, until find a point where trust region method converges and newton doesn't.
    # (if these values stop working, could automate the above process via taking midpoint or centroid.)
    @test_logs (:info, r"converged[\w\s]*trust region") match_mode = :any run_pf([
        35.407000101,
        23.5000028166,
        50.0,
    ])
    # make sure we hit each branch of the trust region method. (code coverage reasons.)
    # TODO hit the "Cauchy step" branch, perhaps by setting a super small starting value of delta.
    @test_logs (:debug, "NR step") (:debug, "Dogleg step") match_mode = :any min_level =
        Logging.Debug run_pf([35.407000101, 23.5000028166, 50.0])
    @test_logs (:debug, "failure") (:debug, "success") match_mode = :any min_level =
        Logging.Debug run_pf([35.407000101, 23.5000028166, 50.0])
    # TODO better tests? i.e. more granularly compare behavior to expected, not just check end result.
    # could check behavior of delta, ie that delta is increased/decreased properly.
end

# TODO: can I create an input such that newton doesn't converge, but iterative_refinement does?
# need jacobian to be ill-conditioned...
