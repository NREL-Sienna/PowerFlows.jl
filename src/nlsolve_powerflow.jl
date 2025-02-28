"""Abstract supertype for all variations of Newton-Raphson."""
abstract type AbstractNewtonRaphsonMethod end
"""Basic, no-frills Newton-Raphson."""
struct SimpleMethod <: AbstractNewtonRaphsonMethod end
"""Use iterative refinement to get better accuracy when solving J(x)*Δx = -f(x)."""
struct RefinementMethod <: AbstractNewtonRaphsonMethod end
"""Trust region method with a dog leg step."""
struct TrustRegionMethod <: AbstractNewtonRaphsonMethod end

"""Cache for non-linear methods
# Fields
-`x::Vector{Float64}`: the current state vector. Used for all methods.
-`r::Vector{Float64}`: the current residual. Used for all methods.
    For `SimpleMethod` and `RefinementMethod`, we solve `J_x Δx = r` in-place,
    so this also stores the step `Δx` at times in those methods.
The remainder of the fields are only used in the `TrustRegionMethod`:
-`r_predict::Vector{Float64}`: the predicted residual at the proposed `x+Δx`,
    under a linear approximation: i.e `J_x⋅(x+p)`.
-`p::Vector{Float64}`: the suggested step `Δx`, selected among `pi`, `p_c`,
    and the dogleg interpolation between the two. The first is chosen when
    `pi` is inside the trust region, the second when both `p_c` and `pi` are outside
    the trust region, and the third when `p_c` is inside and `pi` outside.
    The dogleg step selects the point where the line from `x + p_c` to `x + pi`
    crosses the boundary of the trust region.
-`p_c::Vector{Float64}`: the step to the Cauchy point if the Cauchy point
    lies within the trust region, otherwise a step in that direction.
-`pi::Vector{Float64}`: the step under the Newton-Raphson method."""
struct StateVectorCache
    x::Vector{Float64}
    r::Vector{Float64} # residual
    r_predict::Vector{Float64} # predicted residual
    p::Vector{Float64} # proposed Δx: Cauchy, NR, or dogleg step.
    p_c::Vector{Float64} # Cauchy step
    pi::Vector{Float64} # Newton-Raphson step
end

function StateVectorCache(x0::Vector{Float64}, f0::Vector{Float64})
    x = copy(x0)
    r = copy(f0)
    r_predict = copy(x0)
    p = copy(x0)
    p_c = copy(x0)
    pi = copy(x0)
    return StateVectorCache(x, r, r_predict, p, p_c, pi)
end

"""Reset all entries in a StateVectorCache to `x0` and `f0`, in preparation
for re-using the cache for the next iterative method."""
function resetCache!(
    StateVector::StateVectorCache,
    x0::Vector{Float64},
    f0::Vector{Float64},
)
    copyto!(StateVector.x, x0)
    copyto!(StateVector.r, f0)
    copyto!(StateVector.r_predict, x0)
    copyto!(StateVector.p, x0)
    copyto!(StateVector.p_c, x0)
    copyto!(StateVector.pi, x0)
end

"""Sets `p` equal to the `Δx` by which we should update `x`. Decides
between a Cauchy step (`p = p_c`), Newton-Raphson step (`p = p_i`), and the dogleg
interpolation between the two, based on which fall within the trust region."""
function _dogleg!(p::Vector{Float64},
    p_c::Vector{Float64},
    p_i::Vector{Float64},
    r::Vector{Float64},
    linSolveCache::KLULinSolveCache{Int32},
    Jv::SparseMatrixCSC{Float64, Int32},
    delta::Float64)

    # p_i is newton-raphon step.
    copyto!(p_i, r)
    solve!(linSolveCache, p_i)
    LinearAlgebra.rmul!(p_i, -1.0)

    if norm(p_i) <= delta
        copyto!(p, p_i) # update p: newton-raphson case.
    else
        # using p as a temporary buffer: alias to g for readability
        g = p
        LinearAlgebra.mul!(g, Jv', r)
        p_c .= -norm(g)^2 / norm(Jv * g)^2 .* g # Cauchy point

        if norm(p_c) >= delta
            # p_c outside region => take step of length delta in direction of -g.
            LinearAlgebra.rmul!(g, -delta / norm(g))
            # not needed because g is already an alias for p.
            # copyto!(p, g) # update p: cauchy point case
        else
            # p_c inside region => next point is the spot where the line from p_c to p_i
            # crosses the boundary of the trust region. this is the "dogleg" part.

            # using p_i as temporary buffer: alias to p_diff for readability.
            p_i .-= p_c
            p_diff = p_i

            b = LinearAlgebra.dot(p_c, p_diff)
            a = norm(p_diff)^2
            tau = (-b + sqrt(b^2 - 4a * (norm(p_c)^2 - delta^2))) / (2a)
            p_c .+= tau .* p_diff
            copyto!(p, p_c) # update p: dogleg case.
        end
    end
    return
end

"""Does a single iteration of the `TrustRegionMethod`:
updates the `x` and `r` fields of the `StateVector` and computes
the value of the Jacobian at the new `x`, if needed. Unlike 
`_simple_step`, this has a return value, the updated value of `delta``."""
function _trust_region_step(time_step::Int,
    StateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    Residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    delta::Float64,
    eta::Float64 = DEFAULT_TRUST_REGION_ETA)
    numeric_refactor!(linSolveCache, J.Jv)

    # find proposed next point.
    _dogleg!(
        StateVector.p,
        StateVector.p_c,
        StateVector.pi,
        StateVector.r,
        linSolveCache,
        J.Jv,
        delta,
    )
    StateVector.x .+= StateVector.p

    # use cache.pi as temporary buffer to store old residual
    # to avoid recomputing if we don't change x.
    oldResidual = StateVector.pi
    copyto!(oldResidual, Residual.Rv)
    Residual(StateVector.x, time_step)

    # Ratio of actual to predicted reduction
    LinearAlgebra.mul!(StateVector.r_predict, J.Jv, StateVector.p)
    StateVector.r_predict .+= StateVector.r
    rho =
        (sum(abs2, StateVector.r) - sum(abs2, Residual.Rv)) /
        (sum(abs2, StateVector.r) - sum(abs2, StateVector.r_predict))
    if rho > eta
        # Successful iteration
        StateVector.r .= Residual.Rv
        # we update J here so that if we don't change x (unsuccessful case), we don't re-compute J.
        J(time_step)
    else
        # Unsuccessful: reset x and residual.
        StateVector.x .-= StateVector.p
        copyto!(Residual.Rv, oldResidual)
    end

    # Update size of trust region
    if rho < 0.1 # insufficient improvement
        delta = delta / 2
    elseif rho >= 0.9 # good improvement
        delta = 2 * norm(StateVector.p)
    elseif rho >= 0.5 # so-so improvement
        delta = max(delta, 2 * norm(StateVector.p))
    end
    return delta
end

"""Does a single iteration of either `SimpleMethod` or `RefinementMethod`,
based on the value of `refinement::Bool`. Updates the `r` and `x`
 fields of the `StateVector`, and computes the Jacobian at the new `x`."""
function _simple_step(time_step::Int,
    StateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    Residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    refinement::Bool = false,
    refinement_eps::Float64 = 1e-6)
    # Note: before solve, the R.Rv has the residuals. After solve, R.Rv has the solution (dx). This is
    # because the solve! function modifies the input vector in-place.
    try
        # factorize the numeric object of KLU inplace, while reusing the symbolic object
        numeric_refactor!(linSolveCache, J.Jv)
        # solve for dx in-place
        if !refinement
            solve!(linSolveCache, Residual.Rv)
        else
            Residual.Rv .=
                solve_w_refinement(linSolveCache, J.Jv, Residual.Rv, refinement_eps)
        end
    catch e
        # TODO cook up a test case where Jacobian is singular.
        if e isa LinearAlgebra.SingularException
            @warn("Newton-Raphson hit a point where the Jacobian is singular.")
            fjac2 = similar(J.Jv)
            LinearAlgebra.mul!(fjac2, J.Jv', J.Jv)
            lambda = 1e6 * sqrt(length(StateVector.x) * eps()) * norm(fjac2, 1)
            M = -(fjac2 + lambda * I)
            tempCache = KLULinSolveCache(M) # not reused: just want a minimally-allocating
            # KLU factorization. TODO check if this is faster than Julia's default ldiv.
            full_factor!(tempCache, M)
            solve!(tempCache, Residual.Rv)
        else
            @error("KLU factorization failed: $e")
            # V = _calc_V(data, nlCache.x, time_step)
            # Sbus_result = V .* conj(data.power_network_matrix.data * V)
            return
        end
    end
    # update x
    StateVector.x .-= Residual.Rv
    # update data's fields (the bus angles/voltages) to match x, and update the residual.
    # do this BEFORE updating the Jacobian. The Jacobian computation uses data's fields, not x.
    Residual(StateVector.x, time_step)
    # update jacobian.
    J(time_step)
    return
end

"""Runs the full `SimpleMethod`.
# Keyword arguments:
- `maxIterations::Int`: maximum iterations. Default: $DEFAULT_NR_MAX_ITER.
- `tol::Float64`: tolerance. The iterative search ends when `maximum(abs.(residual)) < tol`.
    Default: $DEFAULT_NR_TOL."""
function _nr_method(time_step::Int,
    StateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    Residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    ::SimpleMethod;
    kwargs...)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    i, converged = 0, false
    while i < maxIterations && !converged
        _simple_step(
            time_step,
            StateVector,
            linSolveCache,
            Residual,
            J,
        )
        converged = norm(Residual.Rv, Inf) < tol
        i += 1
    end
    return converged, i
end

"""Runs the full `RefinementMethod`.
# Keyword arguments:
- `maxIterations::Int`: maximum iterations. Default: $DEFAULT_NR_MAX_ITER.
- `tol::Float64`: tolerance. The iterative search ends when `maximum(abs.(residual)) < tol`.
    Default: $DEFAULT_NR_TOL.
- `refinement_eps::Float64`: run iterative refinement on `J_x Δx = r` until
    `norm(Δx_{i}-Δx_{i+1}, 1)/norm(r,1) < refinement_eps`. Default: 
    $DEFAULT_REFINEMENT_EPS """
function _nr_method(time_step::Int,
    StateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    Residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    ::RefinementMethod;
    kwargs...)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    refinement_eps::Float64 = get(kwargs, :refinement_eps, DEFAULT_REFINEMENT_EPS)
    i, converged = 0, false
    while i < maxIterations && !converged
        _simple_step(
            time_step,
            StateVector,
            linSolveCache,
            Residual,
            J,
            true,
            refinement_eps,
        )
        converged = norm(Residual.Rv, Inf) < tol
        i += 1
    end
    return converged, i
end

"""Runs the full `TrustRegionMethod`.
# Keyword arguments:
- `maxIterations::Int`: maximum iterations. Default: $DEFAULT_NR_MAX_ITER.
- `tol::Float64`: tolerance. The iterative search ends when `maximum(abs.(residual)) < tol`.
    Default: $DEFAULT_NR_TOL.
- `factor::Float64`: the trust region starts out with radius `factor*norm(x_0, 1)`,
    where `x_0` is our initial guess, taken from `data`. Default: $DEFAULT_TRUST_REGION_FACTOR.
- `eta::Float64`: improvement threshold. If the observed improvement in our residual
    exceeds `eta` times the predicted improvement, we accept the new `x_i`.
    Default: $DEFAULT_TRUST_REGION_ETA."""
function _nr_method(time_step::Int,
    StateVector::StateVectorCache,
    linSolveCache::KLULinSolveCache{Int32},
    Residual::ACPowerFlowResidual,
    J::ACPowerFlowJacobian,
    ::TrustRegionMethod;
    kwargs...)
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    factor::Float64 = get(kwargs, :factor, DEFAULT_TRUST_REGION_FACTOR)
    eta::Float64 = get(kwargs, :eta, DEFAULT_TRUST_REGION_ETA)

    delta::Float64 = norm(StateVector.x) > 0 ? factor * norm(StateVector.x) : factor
    i, converged = 0, false
    while i < maxIterations && !converged
        delta = _trust_region_step(
            time_step,
            StateVector,
            linSolveCache,
            Residual,
            J,
            delta,
            eta,
        )
        converged = norm(Residual.Rv, Inf) < tol
        i += 1
    end
    return converged, i
end

function _newton_powerflow(
    pf::ACPowerFlow{NewtonRaphsonACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)
    disable_calc_loss_factors = get(kwargs, :disable_calc_loss_factors, false)
    Residual = ACPowerFlowResidual(data, time_step)
    x0 = calculate_x0(data, time_step)
    Residual(x0, time_step)
    J = PowerFlows.ACPowerFlowJacobian(data, time_step)
    J(time_step)  # we need to fill J with values because at this point it was just initialized

    if sum(Residual.Rv) > 10 * (length(Residual.Rv))
        n_buses = first(size(data.bus_type))
        _, ix = findmax(Residual.Rv)
        bx = ix <= n_buses ? ix : ix - n_buses
        bus_no = data.bus_lookup[bx]
        @warn "Initial guess provided results in a large initial residual. Largest residual at bus $bus_no"
    end

    linSolveCache = KLULinSolveCache(J.Jv)
    symbolic_factor!(linSolveCache, J.Jv)
    StateVector = StateVectorCache(x0, Residual.Rv)

    for method in [SimpleMethod(), RefinementMethod(), TrustRegionMethod()]
        converged, i =
            _nr_method(
                time_step,
                StateVector,
                linSolveCache,
                Residual,
                J,
                method;
                kwargs...,
            )
        if converged
            @info(
                "The NewtonRaphsonACPowerFlow solver converged after $i iterations with method $method"
            )
            V = _calc_V(data, StateVector.x, time_step)
            Sbus_result = V .* conj(data.power_network_matrix.data * V)

            if data.calc_loss_factors && !disable_calc_loss_factors
                calculate_loss_factors(data, J.Jv, time_step)
            end

            return (true, V, Sbus_result)
        end
        @warn("Failed with method $method")
        # reset back to starting point before trying next method.
        # Order here (R then J) is important.
        Residual(x0, time_step)
        J(time_step)
        resetCache!(StateVector, x0, Residual.Rv)
    end

    V = fill(NaN + NaN * im, length(StateVector.x) ÷ 2)
    Sbus_result = fill(NaN + NaN * im, length(StateVector.x) ÷ 2)
    @error(
        "Solver NewtonRaphsonACPowerFlow did not converge with any method."
    )
    return (false, V, Sbus_result)
end

"""
    _calc_V(data::ACPowerFlowData, x::Vector{Float64}, time_step::Int64) -> Vector{Complex{Float64}}

Calculate the results for complex bus voltages from the "x" results of NLSolveACVPowerFlow.
    This is for compatibility with the results of legacy power flow method, ehich returns the vector V instead of the vector x.

# Arguments
- `data::ACPowerFlowData`: The power flow data struct.
- `x::Vector{Float64}`: The results vector from NLSolveACPowerFlow containing voltage magnitudes and angles, as well as active and reactive powers.
- `time_step::Int64`: The time step index for which to calculate the voltages (default is 1).

# Returns
- `Vector{Complex{Float64}}`: A vector of complex bus voltages.

# Details
This function calculates the complex bus voltages based on the bus types:
- REF bus: Voltage magnitude and angle are taken from `data`, because the reference buses maintain the voltage specified by the input data.
- PV bus: Voltage magnitude is taken from `data`, and the angle is taken from `x`. The voltage magnitude is maintained according to the inputs, and the voltage angle is determined in the PF calculation.
- PQ bus: Both voltage magnitude and angle are taken from `x`, as the voltage magnitude and angle are results of the PF calculation for PQ buses.

The state vector `x` is assumed to have 2 values per bus (real and imaginary parts, two of P, Q, Vm (V), Va (θ)).
"""

function _calc_V(
    data::ACPowerFlowData,
    x::Vector{Float64},
    time_step::Int64,
)
    n_buses = length(x) ÷ 2  # Since x has 2 elements per bus (real and imaginary)
    V = zeros(Complex{Float64}, n_buses)
    Vm_data = data.bus_magnitude[:, time_step]
    Va_data = data.bus_angles[:, time_step]
    bus_types = view(data.bus_type, :, time_step)

    # Extract values for Vm and Va from x
    for (ix, bt) in enumerate(bus_types)
        if bt == PSY.ACBusTypes.REF
            # For REF bus, we have active and reactive power
            Vm = Vm_data[ix]
            Va = Va_data[ix]
            V[ix] = Vm * exp(im * Va)
        elseif bt == PSY.ACBusTypes.PV
            # For PV bus, we have reactive power and voltage angle
            Vm = Vm_data[ix]
            Va = x[2 * ix]
            V[ix] = Vm * exp(im * Va)  # Rebuild voltage from magnitude and angle
        elseif bt == PSY.ACBusTypes.PQ
            # For PQ bus, we have voltage magnitude and voltage angle
            Vm = x[2 * ix - 1]
            Va = x[2 * ix]
            V[ix] = Vm * exp(im * Va)  # Rebuild voltage from magnitude and angle
        end
    end

    return V
end

function calculate_x0(data::ACPowerFlowData,
    time_step::Int64)
    n_buses = length(data.bus_type[:, 1])
    x0 = Vector{Float64}(undef, 2 * n_buses)
    state_variable_count = 1
    for (ix, b) in enumerate(data.bus_type[:, time_step])
        if b == PSY.ACBusTypes.REF
            x0[state_variable_count] =
                data.bus_activepower_injection[ix, time_step] -
                data.bus_activepower_withdrawals[ix, time_step]
            x0[state_variable_count + 1] =
                data.bus_reactivepower_injection[ix, time_step] -
                data.bus_reactivepower_withdrawals[ix, time_step]
            state_variable_count += 2
        elseif b == PSY.ACBusTypes.PV
            x0[state_variable_count] =
                data.bus_reactivepower_injection[ix, time_step] -
                data.bus_reactivepower_withdrawals[ix, time_step]
            x0[state_variable_count + 1] = data.bus_angles[ix, time_step]
            state_variable_count += 2
        elseif b == PSY.ACBusTypes.PQ
            x0[state_variable_count] = data.bus_magnitude[ix, time_step]
            x0[state_variable_count + 1] = data.bus_angles[ix, time_step]
            state_variable_count += 2
        else
            throw(ArgumentError("$b not recognized as a bustype"))
        end
    end
    @assert state_variable_count - 1 == n_buses * 2
    return x0
end
