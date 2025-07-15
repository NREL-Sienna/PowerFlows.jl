"""Driver for the LevenbergMaquardtACPowerFlow method: sets up the data 
structures (e.g. residual), runs the powerflow method via calling `_run_powerflow_method` 
on them, then handles post-processing (e.g. loss factors)."""
function _newton_powerflow(
    pf::ACPowerFlow{LevenbergMaquardtACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)
    residual = ACPowerFlowResidual(data, time_step)
    x0 = calculate_x0(data, time_step)
    J = ACPowerFlowJacobian(data, time_step)
    # check if we need to iterate at all: maybe x0 is a solution.
    if _initial_residual!(residual, x0, pf, data, time_step; kwargs...)
        @info("The LevenbergMaquardtACPowerFlow solver converged after 0 iterations.")
        loss_factors_helper!(data, x0, residual, J, time_step)
        return true
    end

    converged, i = _run_powerflow_method(
        time_step,
        x0,
        residual,
        J;
        kwargs...,
    )

    if converged
        @info("The LevenbergMaquardtACPowerFlow solver converged after $i iterations.")
        loss_factors_helper!(data, x0, residual, J, time_step)
        return true
    end
    @error("The LevenbergMaquardtACPowerFlow solver failed to converge.")
    return false
end

# could add DAMPING_INCR and DAMPING_DECR too.
"""Runs the full `LevenbergMaquardtACPowerFlow`.
# Keyword arguments:
- `maxIterations::Int`: maximum iterations. Default: $DEFAULT_NR_MAX_ITER.
- `tol::Float64`: tolerance. The iterative search ends when `norm(abs.(residual)) < tol`.
    Default: $DEFAULT_NR_TOL.
- `λ_0::Float64`: the initial damping parameter. Larger means more damping and a step
    closer to gradient descent. Default: $DEFAULT_λ_0.
- `maxTestλs::Int`: if unable to find a point with a smaller residual vector after
    increasing the damping parameter this many times, end the search with an error.
    Default: $DEFAULT_MAX_TEST_λs"""
function _run_powerflow_method(
    time_step::Int,
    x::Vector{Float64},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian;
    kwargs...,
)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    λ::Float64 = get(kwargs, :λ_0, DEFAULT_λ_0)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    maxTestλs::Int = get(kwargs, :maxTestλs, DEFAULT_MAX_TEST_λs)
    i, converged = 0, false
    residual(x, time_step)
    resSize = dot(residual.Rv, residual.Rv)
    linf = norm(residual.Rv, Inf)
    @debug "initially: sum of squares $(siground(resSize)), L ∞ norm $(siground(linf)), λ = $λ"
    while i < maxIterations && !converged && !isnan(λ)
        λ = update!(x, residual, J, λ, time_step, maxTestλs)
        converged = !isnan(λ) && norm(residual.Rv, Inf) < tol
        i += 1
    end
    if isnan(λ)
        @error "λ is NaN"
    elseif i == maxIterations
        @error "The LevenbergMaquardtACPowerFlow solver didn't coverge in $maxIterations iterations."
    end
    return converged, i
end

# solving (J^T J + λ I) Δx = -J^T r can be numerically unstable
# instead, take the least squares solution to [J; √λ I] * Δx  = [-r; zeros(n)]
# see eg Kaltenbach's master's thesis on levenberg marquardt, p 13.
function betterResidual(
    x::Vector{Float64},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    λ::Float64,
    time_step::Int,
    residualSize::Float64,
)
    residual(x, time_step)
    J(time_step)
    A = vcat(J.Jv, sqrt(λ) * sparse(LinearAlgebra.I, size(J.Jv)))
    b = vcat(-residual.Rv, zeros(size(J.Jv, 1)))
    Δx = A \ b

    temp = residual.Rv .+ J.Jv * Δx

    x_trial = x .+ Δx
    residual(x_trial, time_step)
    newResidualSize = dot(residual.Rv, residual.Rv)

    predicted_reduction = 0.5 * (residualSize - dot(temp, temp) - λ * dot(Δx, Δx))
    actual_reduction = residualSize - newResidualSize
    @assert predicted_reduction > 0
    ρ = actual_reduction / predicted_reduction

    if actual_reduction > 0 && predicted_reduction > 0
        step_size = norm(Δx)
        linf = norm(residual.Rv, Inf)
        @debug "sum of squares $(siground(newResidualSize)), L ∞ norm $(siground(linf)), λ = $(siground(λ)), ||Δx|| = $(siground(step_size))"
        x .+= Δx
        return ρ
    else
        residual(x, time_step)
    end
    return -1.0
end

"""Updates x following to Levenberg-Maquardt, returning the new value of λ.
Current procedure for adjusting λ:
(1) improvement with current λ => λ /= DAMPING_DECR and x += Δx.
(2) else, λ *= DAMPING_INCR.
    (2a) improvement with this λ => keep λ the same and x += Δx.
    (2b) else, go back to (2).
Here, \"improvement with λ\" means: `norm(F(x+Δx), 2) < norm(F(x), 2)`, where
`Δx` is the solution to `(J' * J + λ I) Δx = -J'*F(x)`.
"""
function update!(
    x::Vector{Float64},
    residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    λ::Float64,
    time_step::Int,
    maxTestλs::Int,
)
    residual(x, time_step)
    residualSize = dot(residual.Rv, residual.Rv)
    J(time_step)
    for _ in 1:maxTestλs
        ρ = betterResidual(x, residual, J, λ, time_step, residualSize)
        if ρ > 0
            λ *= max(1 / 3, 1 - (2 * ρ - 1)^3)
            λ = max(λ, 1.0e-8)
            return λ
        end
        λ *= 2.0
    end
    @error "Unable to improve: gave up after increasing damping factor $maxTestλs times."
    return NaN
end
