function _newton_powerflow(pf::ACPowerFlow{<:RobustHomotopyPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)
    Δt_k = get(kwargs, :Δt_k, DEFAULT_Δt_k)
    homHess = HomotopyHessian(data, time_step)
    x = homotopy_x0(data, time_step)
    t_k = 0.0

    # the sparse structure of the Hessian is different at t_k = 0.0 and t_k > 0.0
    # so we need to increase t_k once before we initialize the solver.
    t_k += Δt_k
    homHess(x, t_k, time_step)

    # options: KLUHessianSolver, CholeskyHessianSolver (fastest)
    hSolver = CholeskyHessianSolver(homHess.Hv)
    symbolic_factor!(hSolver, homHess.Hv)

    # Pre-allocate δ buffer to avoid allocations in each Newton step
    δ = zeros(size(x, 1))

    success = true
    while true # go onto next t_k even if search doesn't terminate within max iterations.
        converged_t_k, _ = _second_order_newton(homHess, t_k, time_step, x, hSolver, δ)
        if t_k == 1.0
            success = converged_t_k
            break
        end
        t_k = min(t_k + Δt_k, 1.0)
    end
    if !success
        @warn "RobustHomotopyPowerFlow failed to find a solution"
        return false
    end

    # Note: For homotopy, we don't have a meaningful iteration count at the outer level
    # The finalize function will handle the success case
    return finalize_solver_result!(success, "RobustHomotopyPowerFlow", data, homHess.J, time_step, 0)
end

sig3(x::Float64) = round(x; sigdigits = 3)

function info_helper(homHess::HomotopyHessian, t_k::Float64, F_val::Float64, msg::String)
    r_val = norm(homHess.pfResidual.Rv, Inf)
    @info "t_k = $(sig3(t_k)): $msg, F_k $(sig3(F_val)), residual $(sig3(r_val))"
end

function _second_order_newton(homHess::HomotopyHessian,
    t_k::Float64,
    time_step::Int,
    x::Vector{Float64},
    hSolver::HessianSolver,
    δ::Vector{Float64};
    kwargs...)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)

    i, converged, stop = 0, false, false
    F_val = F_value(homHess, t_k, x, time_step)
    last_tk = t_k == 1.0
    while i < maxIterations && !converged && !stop
        stop = _second_order_newton_step(
            homHess,
            t_k,
            time_step,
            x,
            hSolver,
            δ,
        )
        F_val = F_value(homHess, t_k, x, time_step)
        converged = (last_tk ? check_convergence(homHess.pfResidual, tol) : abs(F_val) < tol)
        i += 1
        if converged
            info_helper(homHess, t_k, F_val, "converged")
        elseif i == maxIterations && !stop
            info_helper(homHess, t_k, F_val, "max iterations")
        end
    end
    return converged, i
end

function _second_order_newton_step(homHess::HomotopyHessian,
    t_k::Float64,
    time_step::Int,
    x::Vector{Float64},
    hSolver::HessianSolver,
    δ::Vector{Float64},
)
    F_val = F_value(homHess, t_k, x, time_step)
    last_step = t_k == 1.0
    homHess(x, t_k, time_step)
    if !last_step && dot(homHess.grad, homHess.grad) < GRAD_ZERO &&
       LinearAlgebra.isposdef(homHess.Hv) # stop case 1: hit local minimum.
        info_helper(homHess, t_k, F_val, "local minimum")
        return true
    end
    modify_and_numeric_factor!(hSolver, homHess.Hv)
    δ .= homHess.grad
    solve!(hSolver, δ)
    δ .*= -1

    # Create objective function
    ϕ = α -> F_value(homHess, t_k, x + α * δ, time_step)

    # Perform line search
    φ_0 = F_val
    dφ_0 = dot(homHess.grad, δ)
    (α_star, F_val) = BackTracking()(ϕ, 1.0, φ_0, dφ_0)

    if !last_step && norm(δ * α_star) < INSUFFICIENT_CHANGE_IN_X
        # stop case 2: slow progress.
        info_helper(homHess, t_k, F_val, "slow progress")
        return true
    end
    x .+= δ * α_star
    return false
end
