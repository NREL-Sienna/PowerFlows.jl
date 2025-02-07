# Make this struct immutable for performance
struct SolverCache
    x::Tx
    xold::Tx
    p::Tx
    g::Tx
end

function SolverCache(df)
    x = copy(df.x_f)
    xold = copy(x)
    p = copy(x)
    g = copy(x)
    return SolverCache(x, xold, p, g)
end

function newton_polar_powerflow!(
    df::OnceDifferentiable,
    initial_x::AbstractArray{T},
    xtol::Real,
    ftol::Real,
    iterations::Integer,
    extended_trace::Bool,
    linesearch,
    linsolve,
    cache=NewtonCache(df)) where {T}
    n = length(initial_x)
    copyto!(cache.x, initial_x)
    value_jacobian!!(df, cache.x)
    check_isfinite(value(df))
    vecvalue = vec(value(df))
    it = 0
    x_converged, f_converged = assess_convergence(initial_x, cache.xold, value(df), NaN, ftol)
    stopped = any(isnan, cache.x) || any(isnan, value(df)) ? true : false

    converged = x_converged || f_converged
    x_ls = copy(cache.x)
    tracing = store_trace || show_trace || extended_trace

    while !stopped && !converged && it < iterations

        it += 1

        if it > 1
            value_jacobian!(df, cache.x)
        end

        try
            linsolve(cache.p, jacobian(df), vec(value(df)))
            rmul!(cache.p, -1)
        catch e
            throw(e)
        end

        copyto!(cache.xold, cache.x)



        if linesearch isa Static
            x_ls .= cache.x .+ cache.p
            value_jacobian!(df, x_ls)
            alpha, ϕalpha = one(real(T)), value(dfo)
        else
            mul!(vec(cache.g), jacobian(df)', vec(value(df)))
            value_gradient!(dfo, cache.x)
            alpha, ϕalpha = linesearch(dfo, cache.x, cache.p, one(real(T)), x_ls, value(dfo), real(dot(cache.g, cache.p)))
        end
        # fvec is here also updated in the linesearch so no need to call f again.
        copyto!(cache.x, x_ls)
        x_converged, f_converged = assess_convergence(cache.x, cache.xold, value(df), xtol, ftol)
        stopped = any(isnan, cache.x) || any(isnan, value(df)) ? true : false

        converged = x_converged || f_converged
    end
    return (converged, x)
end
