
const INSUFFICIENT_CHANGE_IN_X = 10^(-11)
const GRAD_ZERO = 2 * eps()

function mumps_job!(mumps::Mumps, job::Int)
    MUMPS.set_job!(mumps, job)
    MUMPS.invoke_mumps!(mumps)
    return
end

function _newton_powerflow(pf::ACPowerFlow{<:RobustHomotopyPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)

    # TODO this probably isn't the right spot for this: only need to call once,
    # at very start, e.g. when package is loaded. also seems to be a bit picky: I'm
    # frequently getting errors upon exiting the REPL.
    MPI.Init()

    Δt_k = get(kwargs, :Δt_k, DEFAULT_Δt_k)
    homHess = HomotopyHessian(data, time_step)
    t_k_ref = homHess.t_k_ref
    x = homotopy_x0(data, time_step)

    t_k_ref[] += Δt_k
    homHess(x, time_step)

    icntl = deepcopy(MUMPS.default_icntl)
    icntl[4] = 1 # errors only
    mumps = Mumps{Float64}(MUMPS.mumps_symmetric, icntl, MUMPS.default_cntl32)
    MUMPS.associate_matrix!(mumps, homHess.Hv)
    mumps_job!(mumps, ANALYZE)

    success = true
    while true # go onto next t_k even if search doesn't terminate within max iterations.
        converged_t_k, _ = _second_order_newton(homHess, time_step, x, mumps)
        if t_k_ref[] == 1.0
            success = converged_t_k
            break
        end
        t_k_ref[] = min(t_k_ref[] + Δt_k, 1.0)
    end
    if !success
        @warn "RobustHomotopyPowerFlow failed to find a solution"
    end
    return success
end

sig3(x::Float64) = round(x; sigdigits = 3)

function info_helper(homHess::HomotopyHessian, F_val::Float64, msg::String)
    t_k = homHess.t_k_ref[]
    r_val = norm(homHess.pfResidual.Rv, Inf)
    @info "t_k = $(sig3(t_k)): $msg, F_k $(sig3(F_val)), residual $(sig3(r_val))"
end

function _second_order_newton(homHess::HomotopyHessian,
    time_step::Int,
    x::Vector{Float64},
    mumps::Mumps;
    kwargs...)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)

    i, converged, stop = 0, false, false
    F_val = F_value(homHess, x, time_step)
    last_tk = homHess.t_k_ref[] == 1.0
    while i < maxIterations && !converged && !stop
        stop = _second_order_newton_step(
            homHess,
            time_step,
            x,
            mumps,
        )
        F_val = F_value(homHess, x, time_step)
        converged = (last_tk ? norm(homHess.pfResidual.Rv, Inf) : abs(F_val)) < tol
        i += 1
        if converged
            info_helper(homHess, F_val, "converged")
        elseif i == maxIterations && !stop
            info_helper(homHess, F_val, "max iterations")
        end
    end
    return converged, i
end
function _second_order_newton_step(homHess::HomotopyHessian,
    time_step::Int,
    x::Vector{Float64},
    mumps::Mumps{Float64})
    t_k = homHess.t_k_ref[]
    F_val = F_value(homHess, x, time_step)
    last_step = t_k == 1.0
    homHess(x, time_step)
    if !last_step && dot(homHess.grad, homHess.grad) < GRAD_ZERO &&
       LinearAlgebra.isposdef(homHess.Hv) # stop case 1: hit local minimum.
        info_helper(homHess, F_val, "local minimum")
        return true
    end
    # PERF: pass pointers to mumps, so we don't need to run this each time.
    #       MUMPS stores things in COO format, though: just pass a pointer to 
    #       nzval? That's passing a pointer to a struct internals, though... 
    #       See issue #160 in the MUMPS.jl repo.
    MUMPS.associate_matrix!(mumps, homHess.Hv)
    MUMPS.associate_rhs!(mumps, homHess.grad)
    mumps_job!(mumps, FACTOR)
    mumps_job!(mumps, SOLVE)
    # TODO: what if the hessian is singular? is that okay, as long as RHS is still in image?
    δ = -1 * MUMPS.get_solution(mumps)[:, 1]
    mumps_job!(mumps, FACTOR_CLEANUP)

    err = homHess.Hv * δ + homHess.grad
    @assert dot(err, err) < 10 * eps()

    α_star = line_search(x, time_step, homHess, δ)
    if !last_step && norm(δ * α_star) < INSUFFICIENT_CHANGE_IN_X
        # stop case 2: slow progress.
        info_helper(homHess, F_val, "slow progress")
        return true
    end
    x .+= δ * α_star
    return false
end
