
const INSUFFICIENT_CHANGE_IN_X = 10^(-11)
const GRAD_ZERO = 2 * eps()
function _newton_powerflow(pf::ACPowerFlow{<:RobustHomotopyPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)

    # this probably isn't the right spot for this: only need to call once,
    # at very start, e.g. when package is loaded.
    MPI.Init()

    Δt_k = get(kwargs, :Δt_k, DEFAULT_Δt_k)

    homHess = HomotopyHessian(data, 1)
    x = calculate_x0(data, time_step)
    for (bus_ix, bt) in enumerate(get_bus_type(data)[:, time_step])
        if bt == PSY.ACBusTypes.PQ
            x[2 * bus_ix - 1] = 1.0
        end
    end

    t_k = Δt_k
    homHess.t_k_ref[] = t_k
    success::Bool = true
    while true
        success, _ = _second_order_newton(homHess, time_step, x) #, mumps)
        if t_k == 1.0
            break
        end
        t_k = min(1.0, t_k + Δt_k)
        homHess.t_k_ref[] = t_k
    end
    # MUMPS.finalize(mumps)
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
    x::Vector{Float64};
    # mumps::Mumps{Float64};
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
            # mumps
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
    return converged || stop, i
end
function _second_order_newton_step(homHess::HomotopyHessian,
    time_step::Int,
    x::Vector{Float64})
    #mumps::Mumps{Float64})
    t_k = homHess.t_k_ref[]
    F_val = F_value(homHess, x, time_step)
    last_step = t_k == 1.0
    homHess(x, time_step)
    if !last_step && dot(homHess.grad, homHess.grad) < GRAD_ZERO &&
       LinearAlgebra.isposdef(homHess.Hv) # stop case 1: hit local minimum.
        info_helper(homHess, F_val, "local minimum")
        return true
    end

    # TODO: fiugre out how to reuse the mumps object and the symbolic factorization.
    # When I tried naively--see commented out arguments--the first solve was correct,
    # but then future ones were wrong.
    settings = deepcopy(MUMPS.default_icntl)
    settings[4] = 1 # errors only

    mumps = Mumps{Float64}(MUMPS.mumps_symmetric, settings, MUMPS.default_cntl32)
    MUMPS.associate_matrix!(mumps, homHess.Hv)
    MUMPS.factorize!(mumps)
    MUMPS.associate_rhs!(mumps, homHess.grad)
    # TODO: what if the hessian is singular and there's no solution?
    MUMPS.solve!(mumps)
    δ = -1 * MUMPS.get_solution(mumps)[:, 1]

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
