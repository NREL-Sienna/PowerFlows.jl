"""Driver for the LevenbergMarquardtACPowerFlow method: sets up the data 
structures (e.g. residual), runs the powerflow method via calling `_run_powerflow_method` 
on them, then handles post-processing (e.g. loss factors)."""
function _newton_powerflow(
    pf::ACPowerFlow{LevenbergMarquardtACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...,
)
    # setup: common code
    residual = ACPowerFlowResidual(data, time_step)
    x0 = improve_x0(pf, data, residual, time_step)
    if OVERRIDE_x0 && :x0 in keys(kwargs) # for convergence testing
       x0 .= get(kwargs, :x0, x0)
       @warn "Overriding initial guess x0."
       residual(x0, time_step)  # re-calculate residual for new x0: might have changed.
    end
 
    J = ACPowerFlowJacobian(data, time_step)
    J(time_step)
    converged = norm(residual.Rv, Inf) < get(kwargs, :tol, DEFAULT_NR_TOL)
    i = 0

    if !converged
        converged, i = _run_powerflow_method(
            time_step,
            x0,
            residual,
            J;
            kwargs...,
        )
    end

    if converged
        @info("The LevenbergMarquardtACPowerFlow solver converged after $i iterations.")
        if get_calculate_loss_factors(data)
            _calculate_loss_factors(data, J.Jv, time_step)
        end
        if get_calculate_voltage_stability_factors(data)
            _calculate_voltage_stability_factors(data, J.Jv, time_step)
        end
        return true
    end
    @error("The LevenbergMarquardtACPowerFlow solver failed to converge.")
    return false
end

# could add DAMPING_INCR and DAMPING_DECR too.
"""Runs the full `LevenbergMarquardtACPowerFlow`.
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
        λ = update_damping_factor!(x, residual, J, λ, time_step, maxTestλs)
        converged = !isnan(λ) && norm(residual.Rv, Inf) < tol
        i += 1
    end
    if isnan(λ)
        @error "λ is NaN"
    elseif i == maxIterations
        @error "The LevenbergMarquardtACPowerFlow solver didn't coverge in $maxIterations iterations."
    end

    return converged, i
end

# solving (J^T J + λ I) Δx = -J^T r can be numerically unstable
# instead, take the least squares solution to [J; √λ I] * Δx  = [-r; zeros(n)]
# see eg Kaltenbach's master's thesis on levenberg marquardt, p 13.
function compute_error(
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
    b = vcat(-residual.Rv, zeros(size(J.Jv, 2)))
    Δx = A \ b

    temp = residual.Rv .+ J.Jv * Δx

    x_trial = x .+ Δx
    residual(x_trial, time_step)
    newResidualSize = dot(residual.Rv, residual.Rv)

    predicted_reduction = (residualSize - dot(temp, temp) - λ * dot(Δx, Δx))
    actual_reduction = (residualSize - newResidualSize)
    ρ = actual_reduction / predicted_reduction

    if ρ > 1e-4
        x .+= Δx
    end

    return ρ
end


function update_damping_factor!(
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
   
    # This builds a diagonal matrix whose diagonal entries are from Jt*J
    H_diag = LinearAlgebra.Diagonal(sum(abs2, J.Jv; dims=1)[:])

    # Compute eigenvalues of the matrix to use their spacing to update the damping factor
    eigvals = LinearAlgebra.diag(H_diag)
    eigvals_sorted = sort(eigvals, by=x -> abs(x), rev=true)

    # Log-scaled spacings
    log_eigvals = log.(eigvals_sorted)
    spacings = diff(log_eigvals)
    spacing = maximum(spacings)       # Measure differences between consecutive singular values

    # Use the largest spacing as a relaxation factor
    relaxation_factor = maximum(spacing)
    relaxation_factor = max(relaxation_factor, 4.0)  # Prevent very small factors

    test_lambda::Int = 1
    # https://www.sciencedirect.com/science/article/pii/S0142061518336342
    λ = DEFAULT_λ_0*norm(residual.Rv) # this helps with some stability issues
    # delayed gratification seems the more robust approach to the lambda update, and the more standard
    while test_lambda<maxTestλs
        ρ = compute_error(x, residual, J, λ, time_step, residualSize)
        if ρ > 0.75 # good step
            factor = (relaxation_factor>4.0) ? relaxation_factor : 4.0
            λ /= factor
            λ = max(λ, 1.0e-8)
            break
        elseif ρ>=0.25 # okay step
            break
        else # bad step
            factor = (relaxation_factor>0 && relaxation_factor<4.0) ? relaxation_factor : 4.0
            λ *= factor
        end

       test_lambda += 1
    end
    if test_lambda==maxTestλs
        @error "Unable to improve: gave up after increasing damping factor $maxTestλs times."
        λ = NaN
    end

    return λ
end
