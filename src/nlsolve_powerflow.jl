const _NLSOLVE_AC_POWERFLOW_KWARGS =
    Set([:check_reactive_power_limits, :check_connectivity])

# keep around for now for performance comparison reasons.
function _newton_powerflow(
    pf::ACPowerFlow{NLSolveACPowerFlowOld},
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

function _newton_powerflow(
    ::ACPowerFlow{NLSolveACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    # copy-pasted from default options of NLsolve.jl
    xtol::Float64 = 0.0, # NLSolve declares these as real. does the difference matter?
    ftol::Float64 = 1e-8,
    maxIter::Integer = 1_000,
    # unused: added to prevent "no such function" errors from a few tests.
    check_reactive_power_limits = false,
    method = :newton) 

    pf = PolarPowerFlow(data, time_step)
    J = PowerFlows.PolarPowerFlowJacobian(data, pf.x0, time_step)
    x = deepcopy(pf.x0)
    dx = similar(x)

    cache = KLULinSolveCache(J.Jv)
    symbolic_factor!(cache, J.Jv)
    i, converged = 0, false
    while i < maxIter && !converged
        # update jacobian. (when i is 0, already updated by constructor for J.)
        if i > 0
            J(J.Jv, x)
        end

        try
            # factorize the numeric object of KLU inplace, while reusing the symbolic object
            numeric_refactor!(cache, J.Jv)

            # solve for dx in-place
            copyto!(dx, pf.residual)
            solve!(cache, dx)
        catch e
            @error("KLU factorization failed: $e")
            return (converged, x)
        end

        # update x
        x -= dx
        # update residual.
        pf(x) 
        converged = (norm(dx) <= xtol) | (LinearAlgebra.norm(pf.residual, Inf) < ftol)
        i += 1
    end

    # old stuff.
    if !converged
        V = fill(NaN + NaN * im, length(x) ÷ 2)
        Sbus_result = fill(NaN + NaN * im, length(x) ÷ 2)
        @error(
            "Did not converge."
        )
    else
        V = _calc_V(data, x, time_step)
        Sbus_result = V .* conj(data.power_network_matrix.data * V)
    end
    return (converged, V, Sbus_result)
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
