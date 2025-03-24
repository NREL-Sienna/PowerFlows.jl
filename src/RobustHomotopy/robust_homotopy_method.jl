
const INSUFFICIENT_CHANGE_IN_X = 10^(-11)
function _newton_powerflow(pf::ACPowerFlow{<:RobustHomotopyPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)

    # this probably isn't the right spot for this.
    #MPI.Init()
    #mumps = Mumps{Float64}(MUMPS.mumps_symmetric, MUMPS.default_icntl, MUMPS.default_cntl32)

    Δt_k = get(kwargs, :Δt_k, DEFAULT_Δt_k)

    homHess = HomotopyHessian(data, 1)
    x = calculate_x0(data, time_step)
    for (bus_ix, bt) in enumerate(get_bus_type(data)[:, time_step])
        if bt == PSY.ACBusTypes.PQ
            x[2*bus_ix - 1] = 1.0
        end
    end

    t_k = Δt_k
    homHess.t_k_ref[] = t_k
    converging::Bool = true
    while true
        println("t_k = ", round(t_k; sigdigits = 2))
        converging, _ = _second_order_newton(homHess, time_step, x)
        println("converged? $converging")
        g = gradient_value(homHess, x, time_step)
        if LinearAlgebra.isposdef(homHess.Hv) && dot(g, g) < 2*eps()
            println("local minimum")
        end
        if t_k == 1.0
            break
        end
        t_k = min(1.0, t_k + Δt_k)
        homHess.t_k_ref[] = t_k
    end
    # MUMPS.finalize(mumps)
    return converging
end

function _second_order_newton(homHess::HomotopyHessian,
    time_step::Int,
    x::Vector{Float64};
    kwargs...)

    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)

    i, converged, stop = 0, false, false
    F_val = F_value(homHess, x, time_step)
    while i < maxIterations && !converged && !stop
        x_rounded = round.(x, sigdigits=3)
        println("\t$F_val")
        stop = _second_order_newton_step(
            homHess,
            time_step,
            x,
        )
        F_val = F_value(homHess, x, time_step)
        converged = abs(F_val) < tol
        i += 1
    end
    return converged, i
end

function _second_order_newton_step(homHess::HomotopyHessian,
    time_step::Int,
    x::Vector{Float64})


    homHess(x, time_step)
    Hv_dense = Matrix{Float64}(homHess.Hv)
    δ = -1*LinearAlgebra.pinv(Hv_dense)*homHess.grad
    α_star = line_search(x, time_step, homHess, δ)
    # println(round(norm(δ*α_star); sigdigits = 3))
    stop = norm(δ*α_star) < INSUFFICIENT_CHANGE_IN_X
    x .+= δ*α_star
    if stop
        @warn "insufficient change in x: bailing out"
    end

    #= y = homHess.grad
    if all(abs(dot(y, v)) < 10^(-6) for v in eachcol(LinearAlgebra.nullspace(Hv_dense')))
        println("should solve successfully")
    else
        println("gradient fails to be in the image of the hessian")
    end
    # end debugging code
    MUMPS.associate_matrix!(mumps, homHess.Hv)
    MUMPS.factorize!(mumps)
    MUMPS.associate_rhs!(mumps, homHess.grad)
    MUMPS.solve!(mumps)
    x .-= MUMPS.get_solution(mumps)=#
    return stop
end