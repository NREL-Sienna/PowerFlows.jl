"""Cache for non-linear methods.

# Fields
- `x::Vector{Float64}`: the current state vector.
- `r::Vector{Float64}`: the current residual.
- `Δx_nr::Vector{Float64}`: the step under the Newton-Raphson method.
The remainder of the fields are only used in the `TrustRegionACPowerFlow`:
- `r_predict::Vector{Float64}`: the predicted residual at `x+Δx_proposed`,
    under a linear approximation: i.e `J_x⋅(x+Δx_proposed)`.
- `Δx_proposed::Vector{Float64}`: the suggested step `Δx`, selected among `Δx_nr`, 
    `Δx_cauchy`, and the dogleg interpolation between the two. The first is chosen when
    `x+Δx_nr` is inside the trust region, the second when both `x+Δx_cauchy`
    and `x+Δx_nr` are outside the trust region, and the third when `x+Δx_cauchy`
    is inside and `x+Δx_nr` outside. The dogleg step selects the point where the line
    from `x+Δx_cauchy` to `x+Δx_nr` crosses the boundary of the trust region.
- `Δx_cauchy::Vector{Float64}`: the step to the Cauchy point if the Cauchy point
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
function _solve_Δx_nr!(stateVector::StateVectorCache, cache::KLULinSolveCache{J_INDEX_TYPE})
    copyto!(stateVector.Δx_nr, stateVector.r)
    solve!(cache, stateVector.Δx_nr)
    return
end

"""Check error and do refinement."""
function _do_refinement!(stateVector::StateVectorCache,
    A::SparseMatrixCSC{Float64, J_INDEX_TYPE},
    cache::KLULinSolveCache{J_INDEX_TYPE},
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
    linSolveCache::KLULinSolveCache{J_INDEX_TYPE},
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
function _singular_J_fallback(Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE},
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
    Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE},
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
    linSolveCache::KLULinSolveCache{J_INDEX_TYPE},
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

"""Evaluate the Iwamoto objective g(μ) = ‖(1-μ)f₀ + μ²f₁‖² expanded as
g(μ) = (1-μ)²g₀ + 2μ²(1-μ)g₁ + μ⁴g₂ where g₀=‖f₀‖², g₁=f₀ᵀf₁, g₂=‖f₁‖²."""
@inline function _iwamoto_objective(
    μ::Float64, g0::Float64, g1::Float64, g2::Float64,
)::Float64
    om = 1.0 - μ
    μ2 = μ * μ
    return om * om * g0 + 2.0 * μ2 * om * g1 + μ2 * μ2 * g2
end

"""If μ ∈ [IWAMOTO_MU_MIN, IWAMOTO_MU_MAX] and g(μ) < best_g, return the
improved (μ, g(μ)); otherwise return (best_μ, best_g) unchanged."""
@inline function _try_iwamoto_candidate(
    μ::Float64,
    best_μ::Float64,
    best_g::Float64,
    g0::Float64,
    g1::Float64,
    g2::Float64,
)::Tuple{Float64, Float64}
    if IWAMOTO_MU_MIN <= μ <= IWAMOTO_MU_MAX
        gval = _iwamoto_objective(μ, g0, g1, g2)
        if gval < best_g
            return μ, gval
        end
    end
    return best_μ, best_g
end

"""Compute the optimal Iwamoto step multiplier μ ∈ [IWAMOTO_MU_MIN, IWAMOTO_MU_MAX]
by minimizing g(μ) = (1-μ)²g₀ + 2μ²(1-μ)g₁ + μ⁴g₂.

The stationary points satisfy the cubic g'(μ)/2 = 2g₂μ³ - 3g₁μ² + (g₀+2g₁)μ - g₀ = 0.
All real roots are found analytically via the depressed-cubic trigonometric/Cardano form,
and the global minimizer of g over the domain is returned. O(1), zero-allocation."""
function _iwamoto_multiplier(g0::Float64, g1::Float64, g2::Float64)::Float64
    # Initialize best candidate from domain boundaries.
    best_μ = IWAMOTO_MU_MIN
    best_g = _iwamoto_objective(IWAMOTO_MU_MIN, g0, g1, g2)
    best_μ, best_g = _try_iwamoto_candidate(IWAMOTO_MU_MAX, best_μ, best_g, g0, g1, g2)

    # Cubic coefficients: c₃μ³ + c₂μ² + c₁μ + c₀ = 0
    c3 = 2.0 * g2
    c2 = -3.0 * g1
    c1 = g0 + 2.0 * g1
    c0 = -g0

    if abs(c3) < IWAMOTO_DEGENERACY_TOL
        # Degenerate: solve quadratic c₂μ² + c₁μ + c₀ = 0
        if abs(c2) > IWAMOTO_DEGENERACY_TOL
            disc = c1 * c1 - 4.0 * c2 * c0
            if disc >= 0.0
                sq = sqrt(disc)
                for μ in ((-c1 + sq) / (2.0 * c2), (-c1 - sq) / (2.0 * c2))
                    best_μ, best_g = _try_iwamoto_candidate(μ, best_μ, best_g, g0, g1, g2)
                end
            end
        elseif abs(c1) > IWAMOTO_DEGENERACY_TOL
            best_μ, best_g = _try_iwamoto_candidate(-c0 / c1, best_μ, best_g, g0, g1, g2)
        end
        return best_μ
    end

    # Full cubic — depress to t³ + At + B = 0 via μ = t - p/3
    p = c2 / c3
    q = c1 / c3
    c0_n = c0 / c3
    p3 = p / 3.0
    A = q - p * p3
    B = c0_n - q * p3 + 2.0 * p3^3
    Δ = -4.0 * A^3 - 27.0 * B^2

    if Δ > 0.0
        # Three distinct real roots — trigonometric form (A < 0 guaranteed when Δ > 0).
        s = sqrt(-A / 3.0)
        m = 2.0 * s
        arg = clamp(-B / (2.0 * s * s * s), -1.0, 1.0)
        φ3 = acos(arg) / 3.0
        for k in 0:2
            best_μ, best_g = _try_iwamoto_candidate(
                m * cos(φ3 - 2.0 * π * k / 3.0) - p3,
                best_μ, best_g, g0, g1, g2)
        end
    elseif Δ < 0.0
        # One real root — Cardano's formula.
        sqD = sqrt(max(-Δ / 108.0, 0.0))
        best_μ, best_g = _try_iwamoto_candidate(
            cbrt(-B / 2.0 + sqD) + cbrt(-B / 2.0 - sqD) - p3,
            best_μ, best_g, g0, g1, g2)
    else
        # Δ ≈ 0 — repeated roots.
        if abs(A) < IWAMOTO_DEGENERACY_TOL
            # Triple root at t = 0.
            best_μ, best_g = _try_iwamoto_candidate(-p3, best_μ, best_g, g0, g1, g2)
        else
            # Simple root t₁ = 3B/A and double root t₂ = -3B/(2A).
            for t in (3.0 * B / A, -3.0 * B / (2.0 * A))
                best_μ, best_g = _try_iwamoto_candidate(
                    t - p3, best_μ, best_g, g0, g1, g2)
            end
        end
    end

    return best_μ
end

"""Does a single iteration of `NewtonRaphsonACPowerFlow`. Updates the `r` and `x`
 fields of the `stateVector`, and computes the Jacobian at the new `x`."""
function _simple_step(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{J_INDEX_TYPE},
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

"""Does a single iteration of Newton-Raphson with Iwamoto step control.
Computes the Newton step, takes a full trial step, and checks whether the
residual norm decreased. If not, computes an optimal damping multiplier `μ`
and applies a damped step instead. Updates `r`, `x`, and the Jacobian."""
function _iwamoto_step(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{J_INDEX_TYPE},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    refinement_threshold::Float64 = DEFAULT_REFINEMENT_THRESHOLD,
    refinement_eps::Float64 = DEFAULT_REFINEMENT_EPS,
)
    # Save pre-step residual f into stateVector.r
    copyto!(stateVector.r, residual.Rv)
    # Compute Newton step Δx_nr = -J⁻¹f
    _set_Δx_nr!(
        stateVector,
        J,
        linSolveCache,
        NewtonRaphsonACPowerFlow(),
        refinement_threshold,
        refinement_eps,
    )
    # Take full trial step: x += Δx_nr
    stateVector.x .+= stateVector.Δx_nr
    # Evaluate trial residual b = F(x + Δx)
    residual(stateVector.x, time_step)

    # Compute gram scalars for Iwamoto criterion
    g0 = dot(stateVector.r, stateVector.r)
    g1 = dot(stateVector.r, residual.Rv)
    g2 = dot(residual.Rv, residual.Rv)

    if g2 < g0
        # Full step reduced residual — accept it (μ = 1)
        @debug "Iwamoto: full step accepted (g₂/g₀ = $(g2/g0))"
    else
        # Full step did not reduce residual — compute optimal μ
        μ = _iwamoto_multiplier(g0, g1, g2)
        @debug "Iwamoto: damped step μ = $μ (g₂/g₀ = $(g2/g0))"
        # Undo full step and apply damped step
        stateVector.x .-= stateVector.Δx_nr
        stateVector.x .+= μ .* stateVector.Δx_nr
        # Re-evaluate residual at damped point
        residual(stateVector.x, time_step)
    end
    # Update Jacobian at new x
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
function _run_power_flow_method(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{J_INDEX_TYPE},
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
    iwamoto::Bool = get(kwargs, :iwamoto, false)
    i, converged = 1, false
    bus_types = @view get_bus_type(J.data)[:, time_step]
    while i < maxIterations && !converged
        if iwamoto
            _iwamoto_step(
                time_step,
                stateVector,
                linSolveCache,
                residual,
                J,
                refinement_threshold,
                refinement_eps,
            )
        else
            _simple_step(
                time_step,
                stateVector,
                linSolveCache,
                residual,
                J,
                refinement_threshold,
                refinement_eps,
            )
        end
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
function _run_power_flow_method(time_step::Int,
    stateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{J_INDEX_TYPE},
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

function _newton_power_flow(
    pf::ACPowerFlow{T},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...) where {T <: Union{TrustRegionACPowerFlow, NewtonRaphsonACPowerFlow}}

    # setup: common code
    residual, J, x0 = initialize_power_flow_variables(pf, data, time_step; kwargs...)
    converged = norm(residual.Rv, Inf) < get(kwargs, :tol, DEFAULT_NR_TOL)

    i = 0
    if !converged
        linSolveCache = KLULinSolveCache(J.Jv)
        symbolic_factor!(linSolveCache, J.Jv)
        stateVector = StateVectorCache(x0, residual.Rv)
        converged, i = _run_power_flow_method(
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
