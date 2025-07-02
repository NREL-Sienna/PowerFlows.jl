"""Cache for non-linear methods
# Fields
-`x::Vector{Float64}`: the current state vector.
-`r::Vector{Float64}`: the current residual.
-`Δx_nr::Vector{Float64}`: the step under the Newton-Raphson method.
The remainder of the fields are only used in the `TrustRegionACPowerFlow`:
-`r_predict::Vector{Float64}`: the predicted residual at `x+Δx_proposed`,
    under a linear approximation: i.e `J_x⋅(x+Δx_proposed)`.
-`Δx_proposed::Vector{Float64}`: the suggested step `Δx`, selected among `Δx_nr`, 
    `Δx_cauchy`, and the dogleg interpolation between the two. The first is chosen when
    `x+Δx_nr` is inside the trust region, the second when both `x+Δx_cauchy`
    and `x+Δx_nr` are outside the trust region, and the third when `x+Δx_cauchy`
    is inside and `x+Δx_nr` outside. The dogleg step selects the point where the line
    from `x+Δx_cauchy` to `x+Δx_nr` crosses the boundary of the trust region.
-`Δx_cauchy::Vector{Float64}`: the step to the Cauchy point if the Cauchy point
    lies within the trust region, otherwise a step in that direction."""
struct StateVectorCache
    x::Vector{Float64}
    r::Vector{Float64} # residual
    r_predict::Vector{Float64} # predicted residual
    Δx_proposed::Vector{Float64} # proposed Δx: Cauchy, NR, or dogleg step.
    Δx_cauchy::Vector{Float64} # Cauchy step
    Δx_nr::Vector{Float64} # Newton-Raphson step
end

function StateVectorCache(x0::Vector{Float64}, f0::Vector{Float64})
    x = copy(x0)
    r = copy(f0)
    r_predict = copy(x0)
    Δx_proposed = copy(x0)
    Δx_cauchy = copy(x0)
    Δx_nr = copy(x0)
    return StateVectorCache(x, r, r_predict, Δx_proposed, Δx_cauchy, Δx_nr)
end

"""Returns a stand-in matrix for singular J's."""
function _singular_J_fallback(Jv::SparseMatrixCSC{Float64, Int32},
    x::Vector{Float64})
    fjac2 = Jv' * Jv
    lambda = NR_SINGULAR_SCALING * sqrt(length(x) * eps()) * norm(fjac2, 1)
    return -(fjac2 + lambda * LinearAlgebra.I)
end

"""Solve for the Newton-Raphson step, given the factorization object for `J.Jv` 
(if non-singular) or its stand-in (if singular)."""
function _solve_Δx_nr!(stateVector::StateVectorCache, cache::KLULinSolveCache{Int32})
    copyto!(stateVector.Δx_nr, stateVector.r)
    solve!(cache, stateVector.Δx_nr)
    return
end

"""Check error and do refinement."""
function _do_refinement!(stateVector::StateVectorCache,
    A::SparseMatrixCSC{Float64, Int32},
    cache::KLULinSolveCache{Int32},
    refinement_threshold::Float64,
    refinement_eps::Float64,
)
    # use stateVector.r_predict as temporary buffer.
    δ_temp = stateVector.r_predict
    copyto!(δ_temp, A * stateVector.Δx_nr)
    δ_temp .-= stateVector.r
    delta = norm(δ_temp, 1) / norm(stateVector.r, 1)
    if delta > refinement_threshold
        stateVector.Δx_nr .= solve_w_refinement(cache,
            A,
            stateVector.r,
            refinement_eps)
    end
    return
end

"""Sets the Newton-Raphson step. Usually, this is just `J.Jv \\ stateVector.r`, but
`J.Jv` might be singular."""
function _set_Δx_nr!(stateVector::StateVectorCache,
    J::ACPowerFlowJacobian,
    linSolveCache::KLULinSolveCache{Int32},
    solver::ACPowerFlowSolverType,
    refinement_threshold::Float64,
    refinement_eps::Float64)
    try
        numeric_refactor!(linSolveCache, J.Jv)
    catch e
        if e isa LinearAlgebra.SingularException
            @warn("$solver hit a point where the Jacobian is singular.")
            M = _singular_J_fallback(J.Jv, stateVector.x)
            tempCache = KLULinSolveCache(M)
            # creating a new solver cache to solve Mx = stateVector.r once. Default ldiv!
            # might be faster, but singular J's should be rare.
            full_factor!(tempCache, M)
            _solve_Δx_nr!(stateVector, tempCache)
            _do_refinement!(stateVector, M, tempCache, refinement_threshold, refinement_eps)
            LinearAlgebra.rmul!(stateVector.Δx_nr, -1.0)
        else
            @error("KLU factorization failed: $e")
        end
    else
        _solve_Δx_nr!(stateVector, linSolveCache)
        _do_refinement!(
            stateVector,
            J.Jv,
            linSolveCache,
            refinement_threshold,
            refinement_eps,
        )
        LinearAlgebra.rmul!(stateVector.Δx_nr, -1.0)
    end
    return
end

"""Sets `Δx_proposed` equal to the `Δx` by which we should update `x`. Decides
between the Cauchy step `Δx_cauchy`, Newton-Raphson step `Δx_nr`, and the dogleg
interpolation between the two, based on which fall within the trust region."""
function _dogleg!(Δx_proposed::Vector{Float64},
    Δx_cauchy::Vector{Float64},
    Δx_nr::Vector{Float64},
    r::Vector{Float64},
    Jv::SparseMatrixCSC{Float64, Int32},
    delta::Float64,
)
    if norm(Δx_nr) <= delta
        copyto!(Δx_proposed, Δx_nr) # update Δx_proposed: newton-raphson case.
    else
        # using Δx_proposed as a temporary buffer: alias to g for readability
        g = Δx_proposed
        LinearAlgebra.mul!(g, Jv', r)
        Δx_cauchy .= -norm(g)^2 / norm(Jv * g)^2 .* g # Cauchy point

        if norm(Δx_cauchy) >= delta
            # Δx_cauchy outside region => take step of length delta in direction of -g.
            LinearAlgebra.rmul!(g, -delta / norm(g))
            # not needed because g is already an alias for Δx_proposed.
            # copyto!(Δx_proposed, g) # update Δx_proposed: cauchy point case
        else
            # Δx_cauchy inside region => next point is the spot where the line from 
            # Δx_cauchy to Δx_nr crosses the boundary of the trust region.
            # this is the "dogleg" part.

            # using Δx_nr as temporary buffer: alias to Δx_diff for readability.
            Δx_nr .-= Δx_cauchy
            Δx_diff = Δx_nr

            b = LinearAlgebra.dot(Δx_cauchy, Δx_diff)
            a = norm(Δx_diff)^2
            tau = (-b + sqrt(b^2 - 4a * (norm(Δx_cauchy)^2 - delta^2))) / (2a)
            Δx_cauchy .+= tau .* Δx_diff
            copyto!(Δx_proposed, Δx_cauchy) # update Δx_proposed: dogleg case.
        end
    end
    return
end

"""Does a single iteration of the `TrustRegionNRMethod`:
updates the `x` and `r` fields of the `stateVector` and computes
the value of the Jacobian at the new `x`, if needed. Unlike 
`_simple_step`, this has a return value, the updated value of `delta``."""
function _trust_region_step(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    delta::Float64,
    eta::Float64 = DEFAULT_TRUST_REGION_ETA,
)
    _set_Δx_nr!(
        stateVector,
        J,
        linSolveCache,
        TrustRegionACPowerFlow(),
        DEFAULT_REFINEMENT_THRESHOLD,
        DEFAULT_REFINEMENT_EPS,
    )
    _dogleg!(
        stateVector.Δx_proposed,
        stateVector.Δx_cauchy,
        stateVector.Δx_nr,
        stateVector.r,
        J.Jv,
        delta,
    )
    # find proposed next point.
    stateVector.x .+= stateVector.Δx_proposed

    # use cache.Δx_nr as temporary buffer to store old residual
    # to avoid recomputing if we don't change x.
    oldResidual = stateVector.Δx_nr
    copyto!(oldResidual, residual.Rv)
    residual(stateVector.x, time_step)

    # Ratio of actual to predicted reduction
    LinearAlgebra.mul!(stateVector.r_predict, J.Jv, stateVector.Δx_proposed)
    stateVector.r_predict .+= stateVector.r
    rho =
        (sum(abs2, stateVector.r) - sum(abs2, residual.Rv)) /
        (sum(abs2, stateVector.r) - sum(abs2, stateVector.r_predict))
    if rho > eta
        # Successful iteration
        stateVector.r .= residual.Rv
        # we update J here so that if we don't change x (unsuccessful case), we don't re-compute J.
        J(time_step)
    else
        # Unsuccessful: reset x and residual.
        stateVector.x .-= stateVector.Δx_proposed
        copyto!(residual.Rv, oldResidual)
    end

    # Update size of trust region
    if rho < HALVE_TRUST_REGION # rho < 0.1: insufficient improvement
        delta = delta / 2
    elseif rho >= DOUBLE_TRUST_REGION # rho >= 0.9: good improvement
        delta = 2 * norm(stateVector.Δx_proposed)
    elseif rho >= MAX_DOUBLE_TRUST_REGION # rho >= 0.5: so-so improvement
        delta = max(delta, 2 * norm(stateVector.Δx_proposed))
    end
    return delta
end

"""Does a single iteration of `NewtonRaphsonACPowerFlow`. Updates the `r` and `x`
 fields of the `stateVector`, and computes the Jacobian at the new `x`."""
function _simple_step(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    refinement_threshold::Float64 = DEFAULT_REFINEMENT_THRESHOLD,
    refinement_eps::Float64 = DEFAULT_REFINEMENT_EPS,
)
    copyto!(stateVector.r, residual.Rv)
    _set_Δx_nr!(
        stateVector,
        J,
        linSolveCache,
        NewtonRaphsonACPowerFlow(),
        refinement_threshold,
        refinement_eps,
    )
    # update x
    stateVector.x .+= stateVector.Δx_nr
    # update data's fields (the bus angles/voltages) to match x, and update the residual.
    # do this BEFORE updating the Jacobian. The Jacobian computation uses data's fields, not x.
    residual(stateVector.x, time_step)
    # update jacobian.
    J(time_step)
    return
end

"""Runs the full `NewtonRaphsonACPowerFlow`.
# Keyword arguments:
- `maxIterations::Int`: maximum iterations. Default: $DEFAULT_NR_MAX_ITER.
- `tol::Float64`: tolerance. The iterative search ends when `norm(abs.(residual)) < tol`.
    Default: $DEFAULT_NR_TOL.
- `refinement_threshold::Float64`: If the solution to `J_x Δx = r` satisfies
    `norm(J_x Δx - r, 1)/norm(r, 1) > refinement_threshold`, do iterative refinement to
    improve the accuracy. Default: $DEFAULT_REFINEMENT_THRESHOLD.
- `refinement_eps::Float64`: run iterative refinement on `J_x Δx = r` until
    `norm(Δx_{i}-Δx_{i+1}, 1)/norm(r,1) < refinement_eps`. Default: 
    $DEFAULT_REFINEMENT_EPS """
function _run_powerflow_method(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    ::Type{NewtonRaphsonACPowerFlow};
    kwargs...)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    refinement_threshold::Float64 = get(kwargs,
        :refinement_eps,
        DEFAULT_REFINEMENT_THRESHOLD)
    refinement_eps::Float64 = get(kwargs, :refinement_eps, DEFAULT_REFINEMENT_EPS)
    i, converged = 0, false
    while i < maxIterations && !converged
        _simple_step(
            time_step,
            stateVector,
            linSolveCache,
            residual,
            J,
            refinement_threshold,
            refinement_eps,
        )
        converged = norm(residual.Rv, Inf) < tol
        i += 1
    end
    return converged, i
end

"""Runs the full `TrustRegionNRMethod`.
# Keyword arguments:
- `maxIterations::Int`: maximum iterations. Default: $DEFAULT_NR_MAX_ITER.
- `tol::Float64`: tolerance. The iterative search ends when `maximum(abs.(residual)) < tol`.
    Default: $DEFAULT_NR_TOL.
- `factor::Float64`: the trust region starts out with radius `factor*norm(x_0, 1)`,
    where `x_0` is our initial guess, taken from `data`. Default: $DEFAULT_TRUST_REGION_FACTOR.
- `eta::Float64`: improvement threshold. If the observed improvement in our residual
    exceeds `eta` times the predicted improvement, we accept the new `x_i`.
    Default: $DEFAULT_TRUST_REGION_ETA."""
function _run_powerflow_method(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    ::Type{TrustRegionACPowerFlow};
    kwargs...)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    factor::Float64 = get(kwargs, :factor, DEFAULT_TRUST_REGION_FACTOR)
    eta::Float64 = get(kwargs, :eta, DEFAULT_TRUST_REGION_ETA)

    delta::Float64 = norm(stateVector.x) > 0 ? factor * norm(stateVector.x) : factor
    i, converged = 0, false
    while i < maxIterations && !converged
        delta = _trust_region_step(
            time_step,
            stateVector,
            linSolveCache,
            residual,
            J,
            delta,
            eta,
        )
        converged = norm(residual.Rv, Inf) < tol
        i += 1
    end
    return converged, i
end

function _newton_powerflow(
    pf::ACPowerFlow{T},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...) where {T <: Union{TrustRegionACPowerFlow, NewtonRaphsonACPowerFlow}}
    residual = ACPowerFlowResidual(data, time_step)
    x0 = calculate_x0(data, time_step)
    residual(x0, time_step)
    if norm(residual.Rv, 1) > LARGE_RESIDUAL * length(residual.Rv) &&
       get_robust_power_flow(pf)
        improve_x0!(x0, data, time_step, residual)
    else
        @debug "skipping DC powerflow fallback"
    end
    J = PowerFlows.ACPowerFlowJacobian(data, time_step)
    J(time_step)  # we need to fill J with values because at this point it was just initialized

    if sum(abs, residual.Rv) > LARGE_RESIDUAL * length(residual.Rv)
        lg_res, ix = findmax(residual.Rv)
        lg_res_rounded = round(lg_res; sigdigits = 3)
        pow_type = ix % 2 == 1 ? "active" : "reactive"
        bus_ix = div(ix + 1, 2)
        bus_no = axes(data.power_network_matrix, 1)[bus_ix]
        @warn "Initial guess provided results in a large initial residual of $lg_res_rounded. " *
              "Largest residual at bus $bus_no ($bus_ix by matrix indexing; $pow_type power)"
    end

    linSolveCache = KLULinSolveCache(J.Jv)
    symbolic_factor!(linSolveCache, J.Jv)
    stateVector = StateVectorCache(x0, residual.Rv)

    converged, i = _run_powerflow_method(
        time_step,
        stateVector,
        linSolveCache,
        residual,
        J,
        T;
        kwargs...,
    )
    if converged
        @info("The $T solver converged after $i iterations.")
        if data.calculate_loss_factors
            calculate_loss_factors(data, J.Jv, time_step)
        end
        return true
    end
    @error("The $T solver failed to converge.")
    return false
end

"""If initial residual is large, run a DC power flow and see if that gives
a better starting point for angles. Return the original or the result of the DC powerflow,
whichever gives the smaller residual."""
function improve_x0!(x0::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    residual::ACPowerFlowResidual,
)
    @debug "Trying to improve x0 via DC powerflow fallback"
    residualSize = norm(residual.Rv, 1)
    _dc_powerflow_fallback!(data, time_step)
    # is the new starting point better?
    newx0 = calculate_x0(data, time_step)
    residual(newx0, time_step)
    newResidualSize = norm(residual.Rv, 1)
    if newResidualSize < residualSize
        @info "success: DC powerflow fallback yields better x0"
        copyto!(x0, newx0)
        residual(x0, time_step) # re-calculate for new x0.
    else
        @debug "no improvement from DC powerflow fallback"
    end
    return nothing
end

"""Calculate x0 from data."""
function calculate_x0(data::ACPowerFlowData,
    time_step::Int64)
    n_buses = length(data.bus_type[:, 1])
    x0 = Vector{Float64}(undef, 2 * n_buses)
    update_state!(x0, data, time_step)
    # TODO: make contingent on residual size
    enhanced_flat_start!(x0, data, time_step, n_buses)
    return x0
end

function enhanced_flat_start!(
    x0::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    n_buses::Int,
)
    # TODO: use the subnetworks already computed earlier
    subnetworks =
        PNM.find_subnetworks(data.power_network_matrix.adjacency_data, collect(1:n_buses))
    subnetworks = PNM.assign_reference_buses!(
        subnetworks,
        findall(x -> x == PSY.ACBusTypes.REF, data.bus_type[:, time_step]),
    )

    for (ref_bus, subnetwork) in subnetworks
        pv_buses =
            [ix for ix in subnetwork if data.bus_type[ix, time_step] == PSY.ACBusTypes.PV]
        pq_buses =
            [ix for ix in subnetwork if data.bus_type[ix, time_step] == PSY.ACBusTypes.PQ]
        ref_bus_angle = data.bus_angles[ref_bus, time_step]
        (ref_bus_angle == 0.0) || (x0[2 .* [pv_buses; pq_buses]] .= ref_bus_angle)
        (length(pv_buses) == 0) || (length(pq_buses) == 0) ||
            (
                x0[2 .* pq_buses .- 1] .=
                    sum(data.bus_magnitude[pv_buses, time_step]) / length(pv_buses)
            )
    end
    return
end

"""When solving AC power flows, if the initial guess has large residual, we run a DC power 
flow as a fallback. This runs a DC powerflow on `data::ACPowerFlowData` for the given
`time_step`, and writes the solution to `data.bus_angles`."""
function _dc_powerflow_fallback!(data::ACPowerFlowData, time_step::Int)
    # dev note: for DC, we can efficiently solve for all timesteps at once, and we want branch
    # flows. For AC fallback, we're only interested in the current timestep, and no branch flows
    # PERF: if multi-period and multiple time steps have bad initial guesses,
    #       we're re-creating this factorization for each time step. Store it inside
    #       data.aux_network_matrix instead.
    ABA_matrix = data.aux_network_matrix.data
    solver_cache = KLULinSolveCache(ABA_matrix)
    full_factor!(solver_cache, ABA_matrix)
    p_inj =
        data.bus_activepower_injection[data.valid_ix, time_step] -
        data.bus_activepower_withdrawals[data.valid_ix, time_step]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, time_step] .= p_inj
end
