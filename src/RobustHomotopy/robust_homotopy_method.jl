using Optim


function _newton_powerflow(pf::ACPowerFlow{RobustHomoptyPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)

    println("doing a robust homotopy powerflow")

    Δt_k = get(kwargs, :Δt_k, DEFAULT_Δt_k)

    x0 = calculate_x0(data, time_step)
    homFunc = HomotopyFunction(data)
    homGrad = HomotopyGradient(data)
    homHess = HomotopyHessian(data, 1)
    
    # empty constraints.
    lx = Float64[]; ux = Float64[]
    dfc = TwiceDifferentiableConstraints(lx, ux)

    t_k = 0.0
    converging::Bool = true
    while converging
        homGrad(x0)
        homHess(x0)
        # attempt #1: this gives a stack overflow error.
        df = TwiceDifferentiable(homFunc, homGrad, homHess, x0, 0.0, homGrad.grad, homHess.Hv)
        res = optimize(df, dfc, x0, NewtonTrustRegion())

        # attempt #2: this gives a type error about sparse vs non-sparse.
        # res = optimize(homFunc, homGrad, homHess, x0, NewtonTrustRegion())

        # maybe switch to this interface https://jso.dev/NLPModels.jl/dev/api/#Reference-guide
        # and use the trunk solver from https://github.com/JuliaSmoothOptimizers/JSOSolvers.jl

        copyto!(x0, minimizer(res))
        converging = converged(res)
        println("t_k = $t_k: converged $converging to x = $x0")
        clear!(df)
        if t_k == 1.0
            break
        end
        t_k = min(1.0, t_k + Δt_k)
        homFunc.t_k_ref[] = t_k
        homGrad.t_k_ref[] = t_k
        homHess.t_k_ref[] = t_k
    end
    # data should be overwritten with x0
    return converging
end