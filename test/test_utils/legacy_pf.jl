
import PowerFlows

struct LUACPowerFlow <: ACPowerFlowSolverType end  # Only for testing, a basic implementation using LinearAlgebra.lu, allocates a lot of memory

"""This function is to be able to compare the results of the legacy powerflow solver with the new one"""
function _calc_x(
    data::PowerFlows.ACPowerFlowData,
    time_step::Int64,
)
    n_buses = first(size(data.bus_type))
    # mock the expected x format, where the values depend on the type of the bus:
    x = zeros(Float64, 2 * n_buses)
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
        if bt == PSY.ACBusTypes.REF
            # When bustype == REFERENCE PSY.ACBus, state variables are Active and Reactive Power Generated
            x[2 * ix - 1] =
                data.bus_activepower_injection[ix, time_step] -
                data.bus_activepower_withdrawals[ix, time_step]
            x[2 * ix] =
                data.bus_reactivepower_injection[ix, time_step] -
                data.bus_reactivepower_withdrawals[ix, time_step]
        elseif bt == PSY.ACBusTypes.PV
            # When bustype == PV PSY.ACBus, state variables are Reactive Power Generated and Voltage Angle
            x[2 * ix - 1] =
                data.bus_reactivepower_injection[ix, time_step] -
                data.bus_reactivepower_withdrawals[ix, time_step]
            x[2 * ix] = data.bus_angles[ix, time_step]
        elseif bt == PSY.ACBusTypes.PQ
            # When bustype == PQ PSY.ACBus, state variables are Voltage Magnitude and Voltage Angle
            x[2 * ix - 1] = data.bus_magnitude[ix, time_step]
            x[2 * ix] = data.bus_angles[ix, time_step]
        end
    end
    return x
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
    data::PowerFlows.ACPowerFlowData,
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

# this function is for testing purposes only
function _legacy_dSbus_dV(
    V::Vector{Complex{Float64}},
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
)::Tuple{SparseMatrixCSC{Complex{Float64}, Int32}, SparseMatrixCSC{Complex{Float64}, Int32}}
    diagV = SparseArrays.spdiagm(0 => V)
    diagVnorm = SparseArrays.spdiagm(0 => V ./ abs.(V))
    diagIbus = SparseArrays.spdiagm(0 => Ybus * V)
    dSbus_dVm = diagV * conj.(Ybus * diagVnorm) + conj.(diagIbus) * diagVnorm
    dSbus_dVa = 1im * diagV * conj.(diagIbus - Ybus * diagV)
    return dSbus_dVa, dSbus_dVm
end

# this function is for testing purposes only
function _legacy_J(
    dSbus_dVa::SparseMatrixCSC{Complex{Float64}, Int32},
    dSbus_dVm::SparseMatrixCSC{Complex{Float64}, Int32},
    pvpq::Vector{Int64},
    pq::Vector{Int64},
)
    j11 = real(dSbus_dVa[pvpq, pvpq])
    j12 = real(dSbus_dVm[pvpq, pq])
    j21 = imag(dSbus_dVa[pq, pvpq])
    j22 = imag(dSbus_dVm[pq, pq])
    J = sparse([j11 j12; j21 j22])
    return J
end

# legacy NR implementation - here we do not care about allocations, we use this function only for testing purposes
function PowerFlows._newton_powerflow(
    pf::ACPowerFlow{LUACPowerFlow},
    data::PowerFlows.ACPowerFlowData,
    time_step::Int64;
    kwargs...,
)
    # Fetch maxIter and tol from kwargs, or use defaults if not provided
    maxIter = get(kwargs, :maxIter, PowerFlows.DEFAULT_NR_MAX_ITER)
    tol = get(kwargs, :tol, PowerFlows.DEFAULT_NR_TOL)
    i = 0

    Ybus = data.power_network_matrix.data

    # Find indices for each bus type
    ref, pv, pq =
        PowerFlows.bus_type_idx(data, time_step)
    pvpq = [pv; pq]

    npv = length(pv)
    npq = length(pq)
    npvpq = npv + npq

    Vm = data.bus_magnitude[:, time_step]
    Va = data.bus_angles[:, time_step]
    V = zeros(Complex{Float64}, length(Vm))
    V .= Vm .* exp.(1im .* Va)

    # pre-allocate dx
    dx = zeros(Float64, npv + 2 * npq)

    Sbus =
        data.bus_activepower_injection[:, time_step] -
        data.bus_activepower_withdrawals[:, time_step] +
        1im * (
            data.bus_reactivepower_injection[:, time_step] -
            data.bus_reactivepower_withdrawals[:, time_step]
        )

    mis = V .* conj.(Ybus * V) .- Sbus
    F = [real(mis[pvpq]); imag(mis[pq])]

    converged = npvpq == 0

    dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus)
    J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)

    while i < maxIter && !converged
        i += 1

        # using a different factorization than KLU for testing
        factor_J = LinearAlgebra.lu(J)
        dx .= factor_J \ F

        Va[pv] .-= dx[1:npv]
        Va[pq] .-= dx[(npv + 1):(npv + npq)]
        Vm[pq] .-= dx[(npv + npq + 1):(npv + 2 * npq)]
        V .= Vm .* exp.(1im .* Va)
        Vm .= abs.(V)
        Va .= angle.(V)

        mis = V .* conj.(Ybus * V) .- Sbus
        F .= [real(mis[pvpq]); imag(mis[pq])]

        converged = LinearAlgebra.norm(F, Inf) < tol
        if converged
            break
        end

        dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus)
        J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)
    end

    if !converged
        if data.calculate_loss_factors
            data.loss_factors[:, time_step] .= NaN
        end
        @error("The legacy powerflow solver with LU did not converge after $i iterations")
    else
        Sbus_result = V .* conj.(Ybus * V)
        data.bus_magnitude[:, time_step] .= Vm
        data.bus_angles[:, time_step] .= Va
        P_gen = real(Sbus_result) + data.bus_activepower_withdrawals[:, time_step]
        Q_gen = imag(Sbus_result) + data.bus_reactivepower_withdrawals[:, time_step]
        for (ix, bt) in enumerate(data.bus_type[:, time_step])
            if bt == PSY.ACBusTypes.REF
                data.bus_activepower_injection[ix, time_step] = P_gen[ix]
                data.bus_reactivepower_injection[ix, time_step] = Q_gen[ix]
            elseif bt == PSY.ACBusTypes.PV
                data.bus_reactivepower_injection[ix, time_step] = Q_gen[ix]
            end
        end

        if data.calculate_loss_factors
            dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus)
            J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)
            dSbus_dV_ref = collect(real.(hcat(dSbus_dVa[ref, pvpq], dSbus_dVm[ref, pq]))[:])
            J_t = sparse(transpose(J))
            fact = PowerFlows.KLU.klu(J_t)
            lf = fact \ dSbus_dV_ref  # only take the dPref_dP loss factors, ignore dPref_dQ
            data.loss_factors[pvpq, time_step] .= lf[1:npvpq]
            data.loss_factors[ref, time_step] .= 1.0
        end
        if data.calculate_voltage_stability_factors
            LinearAlgebra.__init__()  # to remove warnings
            Gs =
                J[(npvpq + 1):end, (npvpq + 1):end] -
                J[(npvpq + 1):end, 1:npvpq] * inv(collect(J[1:npvpq, 1:npvpq])) *
                J[1:npvpq, (npvpq + 1):end]
            u_1, (σ_1,), v_1, _ = PROPACK.tsvd_irl(Gs; smallest = true, k = 1)
            σ, u, v = PowerFlows.find_sigma_uv(J, npvpq)

            @assert isapprox(σ_1, σ, atol = 1e-6)
            # the sign does not matter
            @assert isapprox(sign(first(u_1)) * u_1, u, atol = 1e-4)
            @assert isapprox(sign(first(v_1)) * v_1, v, atol = 1e-4)

            data.voltage_stability_factors[ref, time_step] .= 0.0
            data.voltage_stability_factors[first(ref), time_step] = σ
            data.voltage_stability_factors[pv, time_step] .= 0.0
            data.voltage_stability_factors[pq, time_step] .= v
        end
        @info("The legacy powerflow solver with LU converged after $i iterations")
    end
    return converged
end
