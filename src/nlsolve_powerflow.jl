struct NLCache{Tx}
    x::Tx
    xold::Tx
    p::Tx
    # g::Tx # only used for more complex linesearch algorithms.
end

function NLCache(x0::Vector{Float64})
    x = copy(x0)
    xold = copy(x)
    p = copy(x)
    # g = copy(x)
    return NLCache(x, xold, p) #, g)
end

function _nr_step(time_step::Int64, nlCache::NLCache,
    linSolveCache::KLULinSolveCache{Int32},
    pf::PolarPowerFlow, J::ACPowerFlowJacobian, strategy::Symbol = :inplace;
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
    J(time_step)
    return
end

function _newton_powerflow(
    pf::ACPowerFlow{HybridACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)

    # Fetch maxIterations, tol from kwargs, or use defaults if not provided
    maxIterations = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol = get(kwargs, :tol, DEFAULT_NR_TOL)

    # when we have multiperiod power flow calculation and want to preserve the ACPowerFlowJacobian instance between time steps:
    J = get(kwargs, :J, PowerFlows.ACPowerFlowJacobian(data, time_step))

    ppf = PolarPowerFlow(data, time_step)
    nlCache = NLCache(ppf.x0)
    linSolveCache = KLULinSolveCache(J.Jv)
    symbolic_factor!(linSolveCache, J.Jv)

    for strategy in [:inplace, :iterative_refinement]
        i, converged = 0, false
        while i < maxIterations && !converged
            _nr_step(time_step, nlCache, linSolveCache, ppf, J, strategy)
            converged = LinearAlgebra.norm(ppf.residual, Inf) < tol
            i += 1
        end
        if converged
            @info(
                "The HybridACPowerFlow solver converged after $i iterations with strategy $strategy"
            )
            V = _calc_V(data, nlCache.x, time_step)
            Sbus_result = V .* conj(data.power_network_matrix.data * V)

            if pf.calc_loss_factors
                num_buses = first(size(data.bus_type))
                ref, pv, pq = bus_type_idx(data, time_step)
                pvpq = vcat(pv, pq)
                pvpq_coords = [
                    x for pair in zip(
                        [2 * x - 1 for x in 1:num_buses if x in pvpq],
                        [2 * x for x in 1:num_buses if x in pvpq],
                    ) for x in pair
                ]
                data.loss_factors[ref, time_step] .= 1.0
                penalty_factors!(
                    J.Jv[pvpq_coords, pvpq_coords],
                    vec(collect(J.Jv[2 .* ref .- 1, pvpq_coords])),
                    view(
                        data.loss_factors,
                        [x for x in 1:num_buses if x in pvpq],
                        time_step,
                    ),
                    [2 * x - 1 for x in 1:length(pvpq)],
                )
            end

            return (true, V, Sbus_result)
        elseif strategy == :inplace
            @warn("Failed with in-place solving. Trying iterative refinement...")
        end
    end
    V = fill(NaN + NaN * im, length(nlCache.x) ÷ 2)
    Sbus_result = fill(NaN + NaN * im, length(nlCache.x) ÷ 2)
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
