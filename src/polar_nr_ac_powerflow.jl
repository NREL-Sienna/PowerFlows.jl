# Heavily based on NLSolve.jl/solvers/newton.jl

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

function assess_convergence(x::Vector{Float64},
    x_previous::Vector{Float64},
    f, # ???
    xtol::Float64,
    ftol::Float64)
    x_converged, f_converged = false, false
    if norm(x - x_previous, Inf) <= xtol
        x_converged = true
    end
    if norm(f, Inf) <= ftol
        f_converged = true
    end

    return x_converged, f_converged
end

function check_isfinite(x::AbstractArray)
    if any(!isfinite, x)
        i = findall(!isfinite, x)
        throw(IsFiniteException(i))
    end
end
#= mutable struct OnceDifferentiable{TF, TDF, TX} <: AbstractObjective
    # function objects (callable)
    f # objective
    df # (partial) derivative of objective
    fdf # objective and (partial) derivative of objective
    # fdf(fx, jx, x) calls j!(jx, x) [updates J], then returns f!(fx, x)

    F::TF # cache for f output
    DF::TDF # cache for df output
    x_f::TX # x used to evaluate f (stored in F)
    x_df::TX # x used to evaluate df (stored in DF)
    f_calls::Vector{Int}
    df_calls::Vector{Int}
end
make_fdf(x, F::Number, f, g!) = (gx, x) -> {j!(gx, x); return f(x)}.
# i.e. update (gx, x) with g!, then return f(x).
=#
# value!!(obj, x): evluates obj at x, stores that value to obj.F, then returns that value
# value(): like value!!, but does not store to obj.F.
function newton_polar_powerflow!(
    obj::OnceDifferentiable,
    initial_x::AbstractArray{T},
    xtol::Real,
    ftol::Real,
    iterations::Integer,
    # extended_trace::Bool,
    linesearch, # gradient descent update step: start with static step length.
    # see https://julianlsolvers.github.io/LineSearches.jl/latest/index.html
    linsolve,
    cache = NewtonCache(obj)) where {T}
    n = length(initial_x)
    copyto!(cache.x, initial_x)
    value_jacobian!!(obj, cache.x) # calls obj.fdf(obj.F, obj.DF, x), then copies x to obj.x_f and obj.x_df
    # but based on context, seems like this is just evaluating J, F, at x0.
    check_isfinite(value(obj)) # value(obj) is getter for obj.F, most recent objective value.
    vecvalue = vec(value(obj))
    it = 0
    x_converged, f_converged =
        assess_convergence(initial_x, cache.xold, value(obj), NaN, ftol)
    stopped = any(isnan, cache.x) || any(isnan, value(obj)) ? true : false

    converged = x_converged || f_converged
    x_ls = copy(cache.x)
    # tracing = store_trace || show_trace || extended_trace

    while !stopped && !converged && it < iterations
        it += 1

        if it > 1
            value_jacobian!(obj, cache.x)
        end

        try
            linsolve(cache.p, jacobian(obj), vec(value(obj))) # jacobian() is getter for DF, most recent Jacobian value.
            rmul!(cache.p, -1)
        catch e
            throw(e)
        end

        copyto!(cache.xold, cache.x)

        # alpha, ϕalpha are unused?
        if linesearch isa Static
            x_ls .= cache.x .+ cache.p
            value_jacobian!(obj, x_ls)
            alpha, ϕalpha = one(real(T)), value(dfo)
        else
            mul!(vec(cache.g), jacobian(obj)', vec(value(obj)))
            value_gradient!(dfo, cache.x)
            alpha, ϕalpha = linesearch(
                dfo,
                cache.x,
                cache.p,
                one(real(T)),
                x_ls,
                value(dfo),
                real(dot(cache.g, cache.p)),
            )
        end
        # fvec is here also updated in the linesearch so no need to call f again.
        copyto!(cache.x, x_ls)
        x_converged, f_converged =
            assess_convergence(cache.x, cache.xold, value(obj), xtol, ftol)
        stopped = any(isnan, cache.x) || any(isnan, value(obj)) ? true : false

        converged = x_converged || f_converged
    end
    return (converged, x)
end
