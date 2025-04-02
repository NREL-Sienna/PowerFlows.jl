# these belong in defintions, but moved here.
# parameters in Wolfe's criteria line search.
const MAX_LINE_SEARCH_ITERS = 20
const MAX_ZOOM_ITERS = 20
const DEFAULT_c_1 = 1e-4
const DEFAULT_c_2 = 0.9
# TODO: chosen randomly. do some research on strategies or fixed values.
const DEFAULT_α_1 = 100.0
const DEFAULT_α_MAX = 1000.0

"""An implementation of Newton-Raphson line search with strong Wolfe conditions."""

function sufficient_decrease(x::Vector{Float64},
    time_step::Int64,
    homInfo::HomotopyHessian,
    δ::Vector{Float64},
    α::Float64,
    c_1::Float64)

    h_x0 = F_value(homInfo, x, time_step)
    h_α = F_value(homInfo, x + α * δ, time_step)

    grad_x0 = gradient_value(homInfo, x, time_step)

    return h_α <= h_x0 + c_1 * dot(grad_x0, δ)
end

function rate_of_change(x_test::Vector{Float64},
    time_step::Int64,
    homInfo::HomotopyHessian,
    δ::Vector{Float64})
    grad = gradient_value(homInfo, x_test, time_step)
    return dot(grad, δ)
end

"""Algorithm 3.6 from Numerical Optimization by Nocedal and Wright. 
Possible point of confusion: α_lo isn't necessarily smaller than α_hi."""
function zoom(x::Vector{Float64},
    time_step::Int64,
    homInfo::HomotopyHessian,
    δ::Vector{Float64},
    α_lo::Float64,
    α_hi::Float64;
    kwargs...)
    c_1 = get(kwargs, :c_1, DEFAULT_c_1)
    c_2 = get(kwargs, :c_2, DEFAULT_c_2)
    h_lo = F_value(homInfo, x + α_lo * δ, time_step)
    j = 0
    while j < MAX_ZOOM_ITERS
        j += 1

        temp_roc = rate_of_change(x + α_lo * δ, time_step, homInfo, δ)
        #println(round.([α_lo, α_hi, temp_roc]; sigdigits = 3))
        # test invariants of the algorithm.
        # (b) α_lo gives sufficient decrease. (Among all step lengths giving 
        # sufficient decrease that we've tried, α_lo gives the smallest value of h.)
        # (c) tangent line at α_lo slopes downward toward α_hi
        @assert temp_roc * (α_hi - α_lo) < 0

        # interpolate. other strategies possible.
        if α_lo != 0.0 && α_hi != 0.0
            α_test = sqrt(α_lo * α_hi)
        else
            α_test = 0.5 * (α_lo + α_hi)
        end
        x_test = x + α_test * δ
        h_test = F_value(homInfo, x_test, time_step)

        if !sufficient_decrease(x, time_step, homInfo, δ, α_test, c_1) || h_test >= h_lo
            # new interval is bottom half.
            α_hi = α_test
        else
            roc_xtest = rate_of_change(x_test, time_step, homInfo, δ)
            if abs(roc_xtest) <= -c_2 * rate_of_change(x, time_step, homInfo, δ)
                # success!
                return α_test
            end
            if roc_xtest * (α_hi - α_lo) >= 0
                # new interval is bottom half, but backward, so condition (c) holds.
                α_lo, α_hi = α_test, α_lo
                h_lo = h_test
            else
                # new interval is top half.
                α_lo = α_test
                h_lo = h_test
            end
        end
    end
    @warn "max zoom iterations exceeded: possibly at a local minimum"
    if α_lo != 0.0 && α_hi != 0.0
        α_test = sqrt(α_lo * α_hi)
    else
        α_test = 0.5 * (α_lo + α_hi)
    end
    return α_test
end

"""Algorithm 3.5 from Numerical Optimization by Nocedal and Wright."""
function line_search(x::Vector{Float64},
    time_step::Int64,
    homInfo::HomotopyHessian,
    δ::Vector{Float64};
    kwargs...)
    α_max = get(kwargs, :α_max, DEFAULT_α_MAX)
    c_1 = get(kwargs, :c_1, DEFAULT_c_1)
    c_2 = get(kwargs, :c_2, DEFAULT_c_2)
    i = 1
    α_i = get(kwargs, :α_max, DEFAULT_α_1)
    α_prev = 0.0
    h_xprev = F_value(homInfo, x, time_step)

    roc_x0 = rate_of_change(x, time_step, homInfo, δ)
    @assert roc_x0 < 0.0
    while i < MAX_LINE_SEARCH_ITERS && α_i < α_max
        x_i = x + α_i * δ
        h_xi = F_value(homInfo, x_i, time_step)
        roc_xi = rate_of_change(x_i, time_step, homInfo, δ)

        cond_i = !sufficient_decrease(x, time_step, homInfo, δ, α_i, c_1)
        cond_ii = h_xi >= h_xprev
        cond_iii = roc_xi >= 0
        @assert cond_i || cond_ii || cond_iii

        if !sufficient_decrease(x, time_step, homInfo, δ, α_i, c_1) ||
           (h_xi >= h_xprev && i > 1)
            return zoom(x, time_step, homInfo, δ, α_prev, α_i, kwargs...)
        end
        if abs(roc_xi) <= -c_2 * roc_x0
            # done, strong Wolfe conditions satisfied.
            return α_i
        end
        if roc_xi >= 0
            return zoom(x, time_step, homInfo, δ, α_i, α_prev, kwargs...)
        end
        α_prev, h_xprev = α_i, h_xi
        # increase α_i. other strategies possible.
        α_i *= 2
        i += 1
    end
    @warn "max line search iterations exceeded"
    return NaN
end
