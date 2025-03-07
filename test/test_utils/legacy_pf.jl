
import PowerFlows

struct LUACPowerFlow <: ACPowerFlowSolverType end  # Only for testing, a basic implementation using LinearAlgebra.lu, allocates a lot of memory

# this function is for testing purposes only
function _legacy_dSbus_dV(
    V::Vector{Complex{Float64}},
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
)::Tuple{SparseMatrixCSC{Complex{Float64}, Int64}, SparseMatrixCSC{Complex{Float64}, Int64}}
    diagV = SparseArrays.spdiagm(0 => V)
    diagVnorm = SparseArrays.spdiagm(0 => V ./ abs.(V))
    diagIbus = SparseArrays.spdiagm(0 => Ybus * V)
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
    data::PowerFlows.ACPowerFlowData,
    time_step::Int64;
    kwargs...,
)
    # Fetch maxIter and tol from kwargs, or use defaults if not provided
    maxIter = get(kwargs, :maxIter, DEFAULT_NR_MAX_ITER)
    tol = get(kwargs, :tol, DEFAULT_NR_TOL)
    i = 0

    Ybus = data.power_network_matrix.data

    # Find indices for each bus type
    pv, pq =
        PowerFlows.bus_type_idx(data, time_step, (PSY.ACBusTypes.PV, PSY.ACBusTypes.PQ))
    pvpq = [pv; pq]

    npv = length(pv)
    npq = length(pq)
    npvpq = npv + npq

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

    dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus)
    J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)

    while i < maxIter && !converged
        i += 1

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
        if converged
            break
        end

        dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus)
        J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)
    end

    if !converged
        V .*= NaN
        Sbus_result = fill(NaN + NaN * im, length(V))
        @error("The legacy powerflow solver with LU did not converge after $i iterations")
    else
        Sbus_result = V .* conj(Ybus * V)
        @info("The legacy powerflow solver with LU converged after $i iterations")
    end
    return (converged, V, Sbus_result)
end
