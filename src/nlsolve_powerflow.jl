using LinearAlgebra
const _NLSOLVE_AC_POWERFLOW_KWARGS =
    Set([:check_reactive_power_limits, :check_connectivity])

# keep around for now for performance comparison reasons.
function _newton_powerflow(
    pf::ACPowerFlow{NLSolveACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;  # not implemented for NLSolve and not used
    nlsolve_kwargs...,
)
    nlsolve_solver_kwargs =
        filter(p -> !(p.first in _NLSOLVE_AC_POWERFLOW_KWARGS), nlsolve_kwargs)

    pf = PolarPowerFlow(data, time_step)
    J = PowerFlows.PolarPowerFlowJacobian(data, pf.x0, time_step)
    df = NLsolve.OnceDifferentiable(pf, J, pf.x0, pf.residual, J.Jv)
    res = NLsolve.nlsolve(df, pf.x0; nlsolve_solver_kwargs...)
    if !res.f_converged
        V = fill(NaN + NaN * im, length(res.zero) ÷ 2)
        Sbus_result = fill(NaN + NaN * im, length(res.zero) ÷ 2)
        @error(
            "The powerflow solver NLSolve did not converge (returned convergence = $(res.f_converged))"
        )
    else
        V = _calc_V(data, res.zero, time_step)
        Sbus_result = V .* conj(data.power_network_matrix.data * V)
    end
    return (res.f_converged, V, Sbus_result)
end

struct NLCache
    x::Vector{Float64}
    xold::Vector{Float64}
    p::Vector{Float64}
    # g::Tx # only used for more complex linesearch algorithms.
end

function NLCache(x0::Vector{Float64})
    x = copy(x0)
    xold = copy(x)
    p = copy(x)
    # g = copy(x)
    return NLCache(x, xold, p) #, g)
end

struct TrustRegionCache
    # Probably could cut down on the number of fields here:
    # I suspect NLSolve only includes some of these for the SolverTrace.
    x::Vector{Float64}
    xold::Vector{Float64}
    r::Vector{Float64} # residual
    r_predict::Vector{Float64} # predicted residual
    p::Vector{Float64} # proposed Δx: Cauchy, NR, or dogleg step.
    p_c::Vector{Float64} # Cauchy step
    pi::Vector{Float64} # Newton-Raphson step
end

function TrustRegionCache(x0::Vector{Float64}, f0::Vector{Float64})
    x = copy(x0)
    xold = copy(x0)
    r = copy(f0)
    r_predict = copy(x0)
    p = copy(x0)
    p_c = copy(x0)
    pi = copy(x0)
    return TrustRegionCache(x, xold, r, r_predict, p, p_c, pi)
end

function dogleg!(p::Vector{Float64}, p_c::Vector{Float64}, p_i::Vector{Float64},
    r::Vector{Float64}, linSolveCache::KLULinSolveCache{Int32}, 
    Jv::SparseMatrixCSC{Float64, Int32}, delta::Float64)
    
    # p_i is newton-raphon step.
    copyto!(p_i, r)
    solve!(linSolveCache, p_i)
    rmul!(p_i, -1.0)

    if norm(p_i) <= delta
        @debug "NR step"
        copyto!(p, p_i) # update p: newton-raphson case.
    else
        # using p as a temporary buffer: alias to g for readability
        g = p
        mul!(g, Jv', r)
        p_c .= -norm(g)^2/norm(Jv*g)^2 .* g # Cauchy point

        if norm(p_c) >= delta
            # p_c outside region => take step of length delta in direction of -g.
            rmul!(g, -delta/norm(g))
            @debug "Cauchy step"
            # not needed because g is already an alias for p.
            # copyto!(p, g) # update p: cauchy point case
        else
            # p_c inside region => next point is the spot where the line from p_c to p_i
            # crosses the boundary of the trust region. this is the "dogleg" part.

            # using p_i as temporary buffer: alias to p_diff for readability.
            @debug "Dogleg step"
            p_i .-= p_c
            p_diff = p_i

            b = dot(p_c, p_diff)
            a = norm(p_diff)^2
            tau = (-b + sqrt(b^2-4a*(norm(p_c)^2 - delta^2)))/(2a)
            p_c .+= tau .* p_diff
            copyto!(p, p_c) # update p: dogleg case.
        end
    end
end

function _trust_region_step(cache::TrustRegionCache, linSolveCache::KLULinSolveCache{Int32},
                pf::PolarPowerFlow, J::PolarPowerFlowJacobian, delta::Float64, eta::Float64 = 1e-4)
    copyto!(cache.xold, cache.x)
    numeric_refactor!(linSolveCache, J.Jv)

    # find proposed next point.
    dogleg!(cache.p, cache.p_c, cache.pi, cache.r, linSolveCache, J.Jv, delta)
    cache.x .+= cache.p

    # TODO this should really be compute-F-but-don't-update-data.
    pf(cache.x)

    # Ratio of actual to predicted reduction
    mul!(cache.r_predict, J.Jv, cache.p)
    cache.r_predict .+= cache.r
    rho = (sum(abs2, cache.r) - sum(abs2, pf.residual)) / (sum(abs2, cache.r) - sum(abs2, cache.r_predict))
    if rho > eta
        # Successful iteration
        @debug "success"
        cache.r .= pf.residual
        J(cache.x) # we update J here: that way, if we don't change x (unsuccessful case), we don't re-compute J.
    else
        # Unsuccessful: reset x and residual.
        @debug "failure"
        cache.x .-= cache.p
        pf(cache.x)
    end

    # Update size of trust region
    if rho < 0.1 # insufficient improvement
        delta = delta/2
    elseif rho >= 0.9 # good improvement
        delta = 2 * norm(cache.p)
    elseif rho >= 0.5 # so-so improvement ??
        delta = max(delta, 2 * norm(cache.p))
    end

end

function _nr_step(nlCache::NLCache, linSolveCache::KLULinSolveCache{Int32},
    pf::PolarPowerFlow, J::PolarPowerFlowJacobian, strategy::Symbol = :inplace;
    refinement_eps::Float64 = 1e-6)
    copyto!(nlCache.xold, nlCache.x)
    try
        # factorize the numeric object of KLU inplace, while reusing the symbolic object
        numeric_refactor!(linSolveCache, J.Jv)
        # solve for dx in-place
        copyto!(nlCache.p, pf.residual)
        if strategy == :inplace
            solve!(linSolveCache, nlCache.p)
        elseif strategy == :iterative_refinement
            nlCache.p .=
                solve_w_refinement(linSolveCache, J.Jv, nlCache.p, refinement_eps)
        end
        rmul!(nlCache.p, -1)
    catch e
        # TODO cook up a test case where Jacobian is singular.
        if e isa SingularException
            @warn("Newton-Raphson hit a point where the Jacobian is singular.")
            fjac2 = J.Jv' * J.Jv
            lambda = 1e6 * sqrt(length(pf.x0) * eps()) * norm(fjac2, 1)
            M = -(fjac2 + lambda * I)
            tempCache = KLULinSolveCache(M) # not reused: just want a minimally-allocating
            # KLU factorization. TODO check if this is faster than Julia's default ldiv.
            full_factor!(tempCache, M)
            copyto!(nlCache.p, pf.residual)
            solve!(tempCache, nlCache.p)
        else
            @error("KLU factorization failed: $e")
            # V = _calc_V(data, nlCache.x, time_step)
            # Sbus_result = V .* conj(data.power_network_matrix.data * V)
            return
        end
    end
    # update x
    nlCache.x .+= nlCache.p
    # update data's fields (the bus angles/voltages) to match x, and update the residual.
    # do this BEFORE updating the Jacobian. The Jacobian computation uses data's fields, not x.
    pf(nlCache.x)
    # update jacobian.
    J(nlCache.x)
    return
end

function calculate_loss_factors(data::ACPowerFlowData, Jv::SparseMatrixCSC{Float64, Int32}, time_step::Int)
    num_buses = first(size(data.bus_type))
    ref, pv, pq = bus_type_idx(data, time_step)
    pvpq = vcat(pv, pq)
    pvpq_coords = [
        x for pair in zip(
            [2 * x - 1 for x in 1:num_buses if x in pvpq],
            [2 * x for x in 1:num_buses if x in pvpq],
        ) for x in pair
    ]
    data.loss_factors[ref, time_step] .= 0.0
    penalty_factors!(
        Jv[pvpq_coords, pvpq_coords],
        vec(collect(Jv[2 .* ref .- 1, pvpq_coords])),
        view(
            data.loss_factors,
            [x for x in 1:num_buses if x in pvpq],
            time_step,
        ),
        [2 * x - 1 for x in 1:length(pvpq)],
    )
end

function _newton_powerflow(
    pf::ACPowerFlow{HybridACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)

    # Fetch maxIterations, tol, and xtol from kwargs, or use defaults if not provided
    maxIterations = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol = get(kwargs, :tol, DEFAULT_NR_TOL)
    xtol = get(kwargs, :xtol, DEFAULT_NR_XTOL)  # default should be 0.0

    ppf = PolarPowerFlow(data, time_step)
    x0 = copy(ppf.x0)
    J = PowerFlows.PolarPowerFlowJacobian(data, ppf.x0, time_step)
    linSolveCache = KLULinSolveCache(J.Jv)
    symbolic_factor!(linSolveCache, J.Jv)

    for strategy in [:inplace, :iterative_refinement]
        nlCache = NLCache(ppf.x0)
        i, converged = 0, false
        while i < maxIterations && !converged
            _nr_step(nlCache, linSolveCache, ppf, J, strategy)
            converged =
                # (norm(nlCache.x - nlCache.xold) <= xtol) |
                (LinearAlgebra.norm(ppf.residual, Inf) < tol)
            i += 1
        end
        if converged
            @info(
                "The HybridACPowerFlow solver converged after $i iterations with strategy $strategy"
            )
            V = _calc_V(data, nlCache.x, time_step)
            Sbus_result = V .* conj(data.power_network_matrix.data * V)

            if pf.calc_loss_factors
                calculate_loss_factors(data, J.Jv, time_step)
            end

            return (true, V, Sbus_result)
        elseif strategy == :inplace
            @warn("Failed with in-place solving. Trying iterative refinement...")
        end
        # reset back to starting point before trying next strategy.
        ppf(x0)
        J(x0)
        nlCache = NLCache(x0)
    end
    
    i2, converged2 = 0, false
    trCache = TrustRegionCache(x0, ppf.residual)
    delta::Float64 = norm(x0) > 0 ? norm(x0) : 1.0
    while i2 < maxIterations && !converged2
        i2 += 1
        @debug "x_$i2: $(trCache.x)"
        @debug "norm(r_$i2): $(norm(ppf.residual))"
        _trust_region_step(trCache, linSolveCache, ppf, J, delta)
        converged2 =
                # (norm(nlCache.x - nlCache.xold) <= xtol) |
                (LinearAlgebra.norm(ppf.residual, Inf) < tol)
    end

    if converged2
        @info(
            "The HybridACPowerFlow solver converged after $i2 iterations with strategy trust region"
        )
        V = _calc_V(data, trCache.x, time_step)
        Sbus_result = V .* conj(data.power_network_matrix.data * V)

        if pf.calc_loss_factors
            calculate_loss_factors(data, J.Jv, time_step)
        end

        return (true, V, Sbus_result)
    end

    V = fill(NaN + NaN * im, length(trCache.x) ÷ 2)
    Sbus_result = fill(NaN + NaN * im, length(trCache.x) ÷ 2)
    @error(
        "Solver HybridACPowerFlow did not converge in $maxIterations iterations with any strategy."
    )
    return (false, V, Sbus_result)
end

"""
    _calc_V(data::ACPowerFlowData, x::Vector{Float64}, time_step::Int64) -> Vector{Complex{Float64}}

Calculate the results for complex bus voltages from the "x" results of NLSolveACVPowerFlow.
    This is for compatibility with the results of KLUACPowerFlow, ehich returns the vector V instead of the vector x.

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
