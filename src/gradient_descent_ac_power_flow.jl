"""
    GradientDescentACPowerFlow <: ACPowerFlowSolverType

An [`ACPowerFlowSolverType`](@ref) that solves AC power flow by minimising ½‖F(x)‖² with
the Adam optimizer and backtracking line search.

# Solver settings (pass via `solver_settings` Dict)
| Key              | Default | Description                          |
|------------------|---------|--------------------------------------|
| `:learning_rate` | `0.01`  | Adam step size η                     |
| `:beta1`         | `0.9`   | 1st moment decay β₁                  |
| `:beta2`         | `0.999` | 2nd moment decay β₂                  |
| `:epsilon`       | `1e-8`  | Numerical stability ε                |

See also: [`ACPowerFlow`](@ref).
"""
struct GradientDescentACPowerFlow <: ACPowerFlowSolverType end

"""
    AdamConfig(; learning_rate=0.01, beta1=0.9, beta2=0.999, epsilon=1e-8)
    AdamConfig(settings::Dict{Symbol, Any})

Configuration for the Adam optimizer used by [`GradientDescentACPowerFlow`](@ref).
"""
Base.@kwdef struct AdamConfig
    learning_rate::Float64 = 0.01
    beta1::Float64 = 0.9
    beta2::Float64 = 0.999
    epsilon::Float64 = 1e-8
end

function AdamConfig(settings::Dict{Symbol, Any})
    return AdamConfig(;
        learning_rate = get(settings, :learning_rate, 0.01),
        beta1 = get(settings, :beta1, 0.9),
        beta2 = get(settings, :beta2, 0.999),
        epsilon = get(settings, :epsilon, 1e-8),
    )
end

"""
    AdamState(n::Int)

Pre-allocated mutable state for the Adam optimizer. All working arrays are allocated once
before the iteration loop to ensure a zero-allocation inner loop.

# Fields
- `g::Vector{Float64}`: gradient vector (length `n`)
- `m::Vector{Float64}`: 1st moment estimate (length `n`)
- `v::Vector{Float64}`: 2nd moment estimate (length `n`)
- `x_save::Vector{Float64}`: saved state vector for line search rollback (length `n`)
- `t::Int`: step counter for bias correction
"""
mutable struct AdamState
    g::Vector{Float64}
    m::Vector{Float64}
    v::Vector{Float64}
    x_save::Vector{Float64}
    t::Int
end

function AdamState(n::Int)
    return AdamState(zeros(n), zeros(n), zeros(n), zeros(n), 0)
end

"""
    compute_gradient!(state::AdamState, J::ACPowerFlowJacobian, R::ACPowerFlowResidual)

Overwrites `state.g` with Jᵀ·F in-place. Uses `mul!` with the `Transpose` wrapper so the
CSC column structure is traversed without materializing a new sparse matrix.
"""
@inline function compute_gradient!(
    state::AdamState,
    J::ACPowerFlowJacobian,
    R::ACPowerFlowResidual,
)
    LinearAlgebra.mul!(state.g, transpose(J.Jv), R.Rv)
    return nothing
end

"""
    adam_step!(x::Vector{Float64}, state::AdamState, cfg::AdamConfig)

Applies one Adam update to `x` in-place. All intermediate quantities are written directly
into `state.m` and `state.v`; bias-corrected values are computed as scalars.
"""
function adam_step!(
    x::Vector{Float64},
    state::AdamState,
    cfg::AdamConfig,
)
    state.t += 1
    t = state.t
    β₁ = cfg.beta1
    β₂ = cfg.beta2
    η = cfg.learning_rate
    ε = cfg.epsilon

    # Bias-correction denominators (scalars — no allocation)
    bc1 = 1.0 - β₁^t
    bc2 = 1.0 - β₂^t

    @inbounds for i in eachindex(x)
        gᵢ = state.g[i]
        state.m[i] = β₁ * state.m[i] + (1.0 - β₁) * gᵢ
        state.v[i] = β₂ * state.v[i] + (1.0 - β₂) * gᵢ * gᵢ
        m̂ᵢ = state.m[i] / bc1
        v̂ᵢ = state.v[i] / bc2
        x[i] -= η * m̂ᵢ / (sqrt(v̂ᵢ) + ε)
    end
    return nothing
end

"""Interpolate `x` between `x_save` (α=0) and current `x` (α=1) at fraction `α`.
Computes `x .= x_save .+ α .* (x .- x_save)` in-place without allocation."""
function _interpolate_x!(x::Vector{Float64}, x_save::Vector{Float64}, α::Float64)
    @inbounds for i in eachindex(x)
        x[i] = x_save[i] + α * (x[i] - x_save[i])
    end
    return nothing
end

const ADAM_BACKTRACK_FACTOR = 0.5
const ADAM_MAX_BACKTRACKS = 10

"""Driver for the GradientDescentACPowerFlow method: sets up the data structures,
runs the Adam-based power flow method with backtracking line search, then handles
post-processing."""
function _newton_power_flow(
    pf::ACPowerFlow{GradientDescentACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...,
)
    residual, J, x0 = initialize_power_flow_variables(pf, data, time_step; kwargs...)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    converged = norm(residual.Rv, Inf) < tol

    i = 0
    if !converged
        cfg = AdamConfig(get_solver_kwargs(pf))
        state = AdamState(length(x0))

        for iter in 1:maxIterations
            # 1. Evaluate residual (already evaluated on first entry)
            if iter > 1
                residual(x0, time_step)
            end
            residual_norm = norm(residual.Rv, Inf)
            if residual_norm < tol
                converged = true
                i = iter
                break
            end

            # 2. Evaluate Jacobian at current x
            J(time_step)

            # 3. Compute gradient in-place: g ← Jᵀ F
            compute_gradient!(state, J, residual)

            # 4. Save current x for line search
            copyto!(state.x_save, x0)
            old_loss = dot(residual.Rv, residual.Rv)

            # 5. Adam parameter update: x ← x − η · Adam(g)
            adam_step!(x0, state, cfg)

            # 6. Write updated x back to ACPowerFlowData fields
            update_data!(data, x0, time_step)

            # 7. Backtracking line search on ½‖F‖²
            residual(x0, time_step)
            new_loss = dot(residual.Rv, residual.Rv)

            α = 1.0
            for _ in 1:ADAM_MAX_BACKTRACKS
                new_loss <= old_loss && break
                α *= ADAM_BACKTRACK_FACTOR
                _interpolate_x!(x0, state.x_save, α)
                update_data!(data, x0, time_step)
                residual(x0, time_step)
                new_loss = dot(residual.Rv, residual.Rv)
            end

            i = iter
        end

        if !converged
            # Check one last time after the final update
            residual(x0, time_step)
            converged = norm(residual.Rv, Inf) < tol
        end
    end

    @info("Final residual size: $(norm(residual.Rv, 2)) L2, $(norm(residual.Rv, Inf)) L∞.")

    if converged
        @info("The GradientDescentACPowerFlow solver converged after $i iterations.")
        if get_calculate_loss_factors(data)
            _calculate_loss_factors(data, J.Jv, time_step)
        end
        if get_calculate_voltage_stability_factors(data)
            _calculate_voltage_stability_factors(data, J.Jv, time_step)
        end
        return true
    end

    @error("The GradientDescentACPowerFlow solver failed to converge after $i iterations. " *
           "Residual L∞: $(norm(residual.Rv, Inf))")
    return false
end
