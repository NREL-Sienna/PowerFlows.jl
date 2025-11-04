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
    d::Vector{Float64}
end

function StateVectorCache(x0::Vector{Float64}, f0::Vector{Float64})
    x = copy(x0)
    r = copy(f0)
    r_predict = copy(x0)
    Δx_proposed = copy(x0)
    Δx_cauchy = copy(x0)
    Δx_nr = copy(x0)
    return StateVectorCache(x, r, r_predict, Δx_proposed, Δx_cauchy, Δx_nr, ones(size(x0)))
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

"""Returns a stand-in matrix for singular J's."""
function _singular_J_fallback(Jv::SparseMatrixCSC{Float64, Int32},
    x::Vector{Float64})
    fjac2 = Jv' * Jv
    lambda = NR_SINGULAR_SCALING * sqrt(length(x) * eps()) * norm(fjac2, 1)
    return -(fjac2 + lambda * LinearAlgebra.I)
end

"""Sets `Δx_proposed` equal to the `Δx` by which we should update `x`. Decides
between the Cauchy step `Δx_cauchy`, Newton-Raphson step `Δx_nr`, and the dogleg
interpolation between the two, based on which fall within the trust region."""
function _dogleg!(Δx_proposed::Vector{Float64},
    Δx_cauchy::Vector{Float64},
    Δx_nr::Vector{Float64},
    r::Vector{Float64},
    Jv::SparseMatrixCSC{Float64, Int32},
    d::Vector{Float64},
    delta::Float64,
)
    nr_norm = wnorm(d, Δx_nr)
    @debug "Trust region: ||Δx_nr|| = $(siground(nr_norm)), δ = $(siground(delta))"

    if nr_norm <= delta
        copyto!(Δx_proposed, Δx_nr) # update Δx_proposed: newton-raphson case.
        @debug "Newton-Raphson step selected (inside trust region)"
    else
        # using Δx_proposed as a temporary buffer: alias to g for readability
        g = Δx_proposed
        LinearAlgebra.mul!(g, Jv', r)
        g .= g ./ d .^ 2
        Δx_cauchy .= -wnorm(d, g)^2 / sum(abs2, Jv * g) .* g # Cauchy point

        cauchy_norm = wnorm(d, Δx_cauchy)
        @debug "Cauchy point: ||Δx_cauchy|| = $(siground(cauchy_norm))"

        if cauchy_norm >= delta
            # Δx_cauchy outside region => take step of length delta in direction of -g.
            LinearAlgebra.rmul!(g, -delta / wnorm(d, g))
            @debug "Cauchy step selected (truncated to trust region boundary)"
            # not needed because g is already an alias for Δx_proposed.
            # copyto!(Δx_proposed, g) # update Δx_proposed: cauchy point case
        else
            # Δx_cauchy inside region => next point is the spot where the line from 
            # Δx_cauchy to Δx_nr crosses the boundary of the trust region.
            # this is the "dogleg" part.

            # using Δx_nr as temporary buffer: alias to Δx_diff for readability.
            Δx_nr .-= Δx_cauchy
            Δx_diff = Δx_nr

            b = wdot(d, Δx_cauchy, d, Δx_diff)
            a = wnorm(d, Δx_diff)^2
            tau = (-b + sqrt(b^2 - 4a * (wnorm(d, Δx_cauchy)^2 - delta^2))) / (2a)
            Δx_cauchy .+= tau .* Δx_diff
            copyto!(Δx_proposed, Δx_cauchy) # update Δx_proposed: dogleg case.
            @debug "Dogleg step selected (τ = $(siground(tau)))"
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
    eta::Float64,
    autoscale::Bool,
)
    old_delta = delta
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
        stateVector.d,
        delta,
    )
    # find proposed next point.
    stateVector.x .+= stateVector.Δx_proposed

    # use cache.Δx_nr as temporary buffer to store old residual
    # to avoid recomputing if we don't change x.
    oldResidual = stateVector.Δx_nr
    copyto!(oldResidual, residual.Rv)
    old_residual_norm = sum(abs2, stateVector.r)
    residual(stateVector.x, time_step)
    new_residual_norm = sum(abs2, residual.Rv)

    # Ratio of actual to predicted reduction
    LinearAlgebra.mul!(stateVector.r_predict, J.Jv, stateVector.Δx_proposed)
    stateVector.r_predict .+= stateVector.r
    rho =
        (sum(abs2, stateVector.r) - sum(abs2, residual.Rv)) /
        (sum(abs2, stateVector.r) - sum(abs2, stateVector.r_predict))

    @debug "Trust region step: ρ = $(siground(rho)), η = $(siground(eta)), ||Δx|| = $(siground(norm(stateVector.Δx_proposed)))"

    if rho > eta
        # Successful iteration
        stateVector.r .= residual.Rv
        residualSize = dot(residual.Rv, residual.Rv)
        linf = norm(residual.Rv, Inf)
        @debug "Step accepted: sum of squares $(siground(residualSize)), L ∞ norm $(siground(linf)), Δ = $(siground(delta)), ||Δx|| = $(siground(norm(stateVector.Δx_proposed)))"
        # we update J here so that if we don't change x (unsuccessful case), we don't re-compute J.
        J(time_step)
        if autoscale
            for i in 1:length(stateVector.x)
                stateVector.d[i] = max(0.1 * stateVector.d[i], norm(view(J.Jv, :, i)))
            end
        else
            @assert all(stateVector.d .== 1.0)
        end
    else
        # Unsuccessful: reset x and residual.
        stateVector.x .-= stateVector.Δx_proposed
        copyto!(residual.Rv, oldResidual)
        @debug "Step rejected: ρ = $(siground(rho)) ≤ η = $(siground(eta))"
    end

    # Update size of trust region
    if rho < HALVE_TRUST_REGION # rho < 0.1: insufficient improvement
        delta = delta / 2
        @debug "Trust region decreased: δ $(siground(old_delta)) → $(siground(delta)) (ρ < $(HALVE_TRUST_REGION))"
    elseif rho >= DOUBLE_TRUST_REGION # rho >= 0.9: good improvement
        delta = 2 * wnorm(stateVector.d, stateVector.Δx_proposed)
        @debug "Trust region increased (good): δ $(siground(old_delta)) → $(siground(delta)) (ρ ≥ $(DOUBLE_TRUST_REGION))"
    elseif rho >= MAX_DOUBLE_TRUST_REGION # rho >= 0.5: so-so improvement
        delta = max(delta, 2 * wnorm(stateVector.d, stateVector.Δx_proposed))
        @debug "Trust region increased (moderate): δ $(siground(old_delta)) → $(siground(delta)) (ρ ≥ $(MAX_DOUBLE_TRUST_REGION))"
    else
        @debug "Trust region unchanged: δ = $(siground(delta))"
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
    validate_vms::Bool = get(
        kwargs,
        :validate_voltages,
        DEFAULT_VALIDATE_VOLTAGES,
    )
    validation_range::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = get(
        kwargs,
        :vm_validation_range,
        DEFAULT_VALIDATION_RANGE,
    )
    i, converged = 1, false
    bus_types = @view get_bus_type(J.data)[:, time_step]
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
        validate_vms && validate_voltages(stateVector.x, bus_types, validation_range, i)
        converged = norm(residual.Rv, Inf) < tol
        if !converged
            i += 1
        end
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
    autoscale::Bool = get(kwargs, :autoscale, DEFAULT_AUTOSCALE)

    if eta > 1.0 || eta < 0.0
        @warn("η = $eta is outside [0, 1]") # eta is set to 2.0 in one test.
    end

    validate_vms::Bool = get(
        kwargs,
        :validate_voltages,
        DEFAULT_VALIDATE_VOLTAGES,
    )
    validation_range::NamedTuple{(:min, :max), Tuple{Float64, Float64}} = get(
        kwargs,
        :vm_validation_range,
        DEFAULT_VALIDATION_RANGE,
    )

    if autoscale
        for i in 1:length(stateVector.x)
            stateVector.d[i] = norm(view(J.Jv, :, i))
            if stateVector.d[i] == 0.0
                stateVector.d[i] = 1.0
            end
        end
    end

    delta::Float64 = norm(stateVector.x) > 0 ? factor * norm(stateVector.x) : factor
    i, converged = 0, false
    residualSize = dot(residual.Rv, residual.Rv)
    linf = norm(residual.Rv, Inf)
    @debug "initially: sum of squares $(siground(residualSize)), L ∞ norm $(siground(linf)), Δ $(siground(delta))"

    bus_types = @view get_bus_type(J.data)[:, time_step]
    while i < maxIterations && !converged
        delta = _trust_region_step(
            time_step,
            stateVector,
            linSolveCache,
            residual,
            J,
            delta,
            eta,
            autoscale,
        )
        validate_vms && validate_voltages(stateVector.x, bus_types, validation_range, i)
        converged = norm(residual.Rv, Inf) < tol
        if !converged
            i += 1
        end
    end
    return converged, i
end

function _newton_powerflow(
    pf::ACPowerFlow{T},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...) where {T <: Union{TrustRegionACPowerFlow, NewtonRaphsonACPowerFlow}}

    # setup: common code
    residual, J, x0 = initialize_powerflow_variables(pf, data, time_step; kwargs...)
    converged = norm(residual.Rv, Inf) < get(kwargs, :tol, DEFAULT_NR_TOL)

    i = 0
    if !converged
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
    end
    @info("Final residual size: $(norm(residual.Rv, 2)) L2, $(norm(residual.Rv, Inf)) L∞.")

    if converged
        @info("The $T solver converged after $i iterations.")
        if get_calculate_loss_factors(data)
            _calculate_loss_factors(data, J.Jv, time_step)
        end
        if get_calculate_voltage_stability_factors(data)
            _calculate_voltage_stability_factors(data, J.Jv, time_step)
        end
        return true
    end
    @error("The $T solver failed to converge.")
    return false
end
