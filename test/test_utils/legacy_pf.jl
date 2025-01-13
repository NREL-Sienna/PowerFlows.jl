
import PowerFlows

struct LUACPowerFlow <: ACPowerFlowSolverType end  # Only for testing, a basic implementation using LinearAlgebra.lu, allocates a lot of memory

# this function is for testing purposes only
function _legacy_dSbus_dV(
    V::Vector{Complex{Float64}},
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
)
    diagV = LinearAlgebra.Diagonal(V)
    diagVnorm = LinearAlgebra.Diagonal(V ./ abs.(V))
    diagIbus = LinearAlgebra.Diagonal(Ybus * V)
    dSbus_dVm = diagV * conj.(Ybus * diagVnorm) + conj.(diagIbus) * diagVnorm
    dSbus_dVa = 1im * diagV * conj.(diagIbus - Ybus * diagV)
    return dSbus_dVa, dSbus_dVm
end

# this function is for testing purposes only
function _legacy_J(
    dSbus_dVa::SparseMatrixCSC{Complex{Float64}, Int64},
    dSbus_dVm::SparseMatrixCSC{Complex{Float64}, Int64},
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
function _newton_powerflow(
    pf::ACPowerFlow{PowerFlows.LUACPowerFlow},
    data::PowerFlows.ACPowerFlowData;
    time_step::Int64 = 1,
    kwargs...,
)
    # Fetch maxIter and tol from kwargs, or use defaults if not provided
    maxIter = get(kwargs, :maxIter, DEFAULT_NR_MAX_ITER)
    tol = get(kwargs, :tol, DEFAULT_NR_TOL)
    i = 0

    Ybus = data.power_network_matrix.data

    # Find indices for each bus type
    #ref = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.REF, data.bus_type)
    pv = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PV, data.bus_type)
    pq = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PQ, data.bus_type)
    pvpq = [pv; pq]

    #nref = length(ref)
    npv = length(pv)
    npq = length(pq)
    npvpq = npv + npq
    n_buses = length(data.bus_type)

    Vm = data.bus_magnitude[:, time_step]
    # prevent unfeasible starting values for Vm; for pv and ref buses we cannot do this:
    Vm[pq] .= clamp.(Vm[pq], 0.9, 1.1)
    Va = data.bus_angles[:, time_step]
    V = zeros(Complex{Float64}, length(Vm))
    V .= Vm .* exp.(1im .* Va)

    # pre-allocate dx
    dx = zeros(Float64, npv + 2 * npq)

    Ybus = data.power_network_matrix.data

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

    while i < maxIter && !converged
        i += 1
        dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus)
        J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)

        # using a different factorization that KLU for testing
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
    end

    if !converged
        @error("The powerflow solver with KLU did not converge after $i iterations")
    else
        @info("The powerflow solver with KLU converged after $i iterations")
    end
    Sbus_result = V .* conj(Ybus * V)
    return (converged, V, Sbus_result)
end
