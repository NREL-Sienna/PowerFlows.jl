
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
    Ybus::SparseMatrixCSC{Complex{Float32}, Int64},
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
    size_J::Union{Nothing, Int64} = nothing,
)
    j11 = real(dSbus_dVa[pvpq, pvpq])
    j12 = real(dSbus_dVm[pvpq, pq])
    j21 = imag(dSbus_dVa[pq, pvpq])
    j22 = imag(dSbus_dVm[pq, pq])
    J = sparse([j11 j12; j21 j22])
    if !isnothing(size_J)
        I, J, V = findnz(J)
        J = sparse(I, J, V, size_J, size_J)
    end
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

    lcc_p_set = data.lcc_P_set[:, time_step]
    lcc_x_t_i = data.x_t_i[:, time_step]
    lcc_x_t_j = data.x_t_j[:, time_step]
    lcc_I_dc_i = data.lcc_P_set[:, time_step] #[0.01 for _ in lcc_p_set] # This is a placeholder, should be set to the actual DC current value
    lcc_I_dc_j = -data.lcc_P_set[:, time_step]  #[-0.01 for _ in lcc_p_set] # This is a placeholder, should be set to the actual DC current value
    lcc_alpha_i = data.alpha_i[:, time_step]
    lcc_alpha_j = data.alpha_j[:, time_step]
    lcc_t_i = data.t_i[:, time_step]
    lcc_t_j = data.t_j[:, time_step]
    lcc_i = data.lcc_i
    lcc_j = data.lcc_j
    R = data.R[:, time_step]

    n_lcc = length(lcc_p_set)

    Ybus_lcc = _ybus_lcc(
        Ybus, Vm, lcc_t_i, lcc_t_j, lcc_alpha_i, lcc_alpha_j,
        lcc_I_dc_i, lcc_I_dc_j, lcc_x_t_i, lcc_x_t_j, lcc_i, lcc_j,
    )

    # pre-allocate dx
    dx = zeros(Float64, npv + 2 * npq + 4 * n_lcc)

    Sbus =
        data.bus_activepower_injection[:, time_step] -
        data.bus_activepower_withdrawals[:, time_step] +
        1im * (
            data.bus_reactivepower_injection[:, time_step] -
            data.bus_reactivepower_withdrawals[:, time_step]
        )

    F_lcc = _f_lcc(R, lcc_t_i, lcc_t_j, lcc_alpha_i, lcc_alpha_j, lcc_I_dc_i, lcc_I_dc_j,
        lcc_x_t_i, lcc_x_t_j, lcc_p_set, Vm, data.lcc_i, data.lcc_j)

    @show F_lcc

    display(Ybus)

    display(Ybus .+ Ybus_lcc)

    I_ac_lcc = Ybus_lcc * V

    I_ac_lcc_i = I_ac.(lcc_t_i, lcc_alpha_i, lcc_I_dc_i, lcc_x_t_i, Vm[lcc_i], Va[lcc_i])
    I_ac_lcc_j = I_ac.(lcc_t_j, lcc_alpha_j, lcc_I_dc_j, lcc_x_t_j, Vm[lcc_j], Va[lcc_j])
    I_ac_lcc[lcc_i] .= I_ac_lcc_i
    I_ac_lcc[lcc_j] .= I_ac_lcc_j
    @show I_ac_lcc_i
    @show I_ac_lcc_j
    @show I_ac_lcc

    mis = V .* conj.((Ybus .+ Ybus_lcc) * V) .- Sbus .- V .* conj.(I_ac_lcc)
    F = [real(mis[pvpq]); imag(mis[pq]); F_lcc]

    size_J = length(F)

    converged = npvpq == 0

    dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus .+ Ybus_lcc)
    J =
        _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq, size_J) .+
        _legacy_J_lcc(pvpq, pq, Vm, lcc_t_i, lcc_t_j, lcc_I_dc_i, lcc_I_dc_j,
            lcc_alpha_i, lcc_alpha_j, lcc_x_t_i, lcc_x_t_j, size_J,
            lcc_i, lcc_j)

    while i < maxIter && !converged
        i += 1

        # using a different factorization than KLU for testing
        factor_J = LinearAlgebra.lu(J)
        dx .= factor_J \ F

        Va[pv] .-= dx[1:npv]
        Va[pq] .-= dx[(npv + 1):(npv + npq)]
        Vm[pq] .-= dx[(npv + npq + 1):(npv + 2 * npq)]
        d_lcc = -dx[(npv + 2 * npq + 1):(npv + 2 * npq + 4 * n_lcc)]
        V .= Vm .* exp.(1im .* Va)
        Vm .= abs.(V)
        Va .= angle.(V)

        inds1 = [4 * (k - 1) + 1 for k in 1:n_lcc]
        inds2 = [4 * (k - 1) + 2 for k in 1:n_lcc]
        inds3 = [4 * (k - 1) + 3 for k in 1:n_lcc]
        inds4 = [4 * (k - 1) + 4 for k in 1:n_lcc]
        lcc_t_i .+= d_lcc[inds1]
        lcc_t_j .+= d_lcc[inds2]
        lcc_alpha_i .+= d_lcc[inds3]
        lcc_alpha_j .+= d_lcc[inds4]

        clamp!(lcc_t_i, 0.5, 1.5)
        clamp!(lcc_t_j, 0.5, 1.5)  # t_j is wrong (very high number)

        # lcc_alpha_i .= 0.0
        # lcc_alpha_j .= π / 4

        I_ac_lcc_i =
            I_ac.(lcc_t_i, lcc_alpha_i, lcc_I_dc_i, lcc_x_t_i, Vm[lcc_i], Va[lcc_i])
        I_ac_lcc_j =
            I_ac.(lcc_t_j, lcc_alpha_j, lcc_I_dc_j, lcc_x_t_j, Vm[lcc_j], Va[lcc_j])
        I_ac_lcc[lcc_i] .= I_ac_lcc_i
        I_ac_lcc[lcc_j] .= I_ac_lcc_j
        @show I_ac_lcc_i
        @show I_ac_lcc_j
        @show I_ac_lcc

        I_ac_ = Ybus_lcc * V
        @show abs.(V)
        @show I_ac_
        # lcc_I_dc_i .= abs.(sqrt(6) / π .* I_ac[lcc_i]) .* sign.(real.(I_ac[lcc_i]))
        # lcc_I_dc_j .= abs.(sqrt(6) / π .* I_ac[lcc_j]) .* sign.(real.(I_ac[lcc_j]))

        @show lcc_I_dc_i, lcc_I_dc_j

        Ybus_lcc = _ybus_lcc(
            Ybus, Vm, lcc_t_i, lcc_t_j, lcc_alpha_i, lcc_alpha_j,
            lcc_I_dc_i, lcc_I_dc_j, lcc_x_t_i, lcc_x_t_j, lcc_i, lcc_j,
        )

        @show lcc_t_i, lcc_t_j, lcc_alpha_i, lcc_alpha_j

        # @show P_CSC(lcc_t_i, lcc_alpha_i, lcc_I_dc, lcc_x_t_i, Vm[data.lcc_i])
        # @show P_CSC(lcc_t_j, lcc_alpha_j, -lcc_I_dc, lcc_x_t_j, Vm[data.lcc_j])

        F_lcc =
            _f_lcc(R, lcc_t_i, lcc_t_j, lcc_alpha_i, lcc_alpha_j, lcc_I_dc_i, lcc_I_dc_j,
                lcc_x_t_i, lcc_x_t_j, lcc_p_set, Vm, data.lcc_i, data.lcc_j)

        @show F_lcc

        mis = V .* conj.((Ybus .+ Ybus_lcc) * V) .- Sbus .- V .* conj.(I_ac_lcc)
        F .= [real(mis[pvpq]); imag(mis[pq]); F_lcc]

        @show real(mis[pvpq]), imag(mis[pq])

        converged = LinearAlgebra.norm(F, Inf) < tol
        if converged
            break
        end

        dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus .+ Ybus_lcc)
        J =
            _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq, size_J) .+
            _legacy_J_lcc(pvpq, pq, Vm, lcc_t_i, lcc_t_j, lcc_I_dc_i, lcc_I_dc_j,
                lcc_alpha_i, lcc_alpha_j, lcc_x_t_i, lcc_x_t_j, size_J,
                lcc_i, lcc_j)
    end

    if !converged
        if data.calculate_loss_factors
            data.loss_factors[:, time_step] .= NaN
        end
        @error("The legacy powerflow solver with LU did not converge after $i iterations")
    else
        Sbus_result = V .* conj.((Ybus .+ Ybus_lcc) * V)
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

        if PowerFlows.get_calculate_loss_factors(data)
            dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, Ybus)
            J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)
            dSbus_dV_ref = collect(real.(hcat(dSbus_dVa[ref, pvpq], dSbus_dVm[ref, pq]))[:])
            J_t = sparse(transpose(J))
            fact = PowerFlows.KLU.klu(J_t)
            lf = fact \ dSbus_dV_ref  # only take the dPref_dP loss factors, ignore dPref_dQ
            data.loss_factors[pvpq, time_step] .= lf[1:npvpq]
            data.loss_factors[ref, time_step] .= 1.0
        end
        if PowerFlows.get_calculate_voltage_stability_factors(data)
            σ, u, v = PowerFlows._singular_value_decomposition(J, npvpq)
            data.voltage_stability_factors[ref, time_step] .= 0.0
            data.voltage_stability_factors[first(ref), time_step] = σ
            data.voltage_stability_factors[pv, time_step] .= 0.0
            data.voltage_stability_factors[pq, time_step] .= v
        end
        @info("The legacy powerflow solver with LU converged after $i iterations")
    end
    return converged
end

function _ybus_lcc(
    Ybus,
    Vm,
    lcc_t_i,
    lcc_t_j,
    lcc_alpha_i,
    lcc_alpha_j,
    lcc_I_dc_i,
    lcc_I_dc_j,
    lcc_x_i,
    lcc_x_j,
    lcc_i,
    lcc_j,
)
    I = Int[]
    J = Int[]
    V = ComplexF32[]

    for k in eachindex(lcc_i)
        i = lcc_i[k]
        j = lcc_j[k]
        α_i = lcc_alpha_i[k]
        α_j = lcc_alpha_j[k]
        x_i = lcc_x_i[k]
        x_j = lcc_x_j[k]
        t_i = lcc_t_i[k]
        t_j = lcc_t_j[k]
        I_dc_i = lcc_I_dc_i[k]
        I_dc_j = lcc_I_dc_j[k]

        Y_i = Y_val(t_i, α_i, I_dc_i, x_i, Vm[i])
        @show t_j, α_j, I_dc_j, x_j, Vm[j]
        Y_j = Y_val(t_j, α_j, I_dc_j, x_j, Vm[j])

        push!(I, i)
        push!(J, i)
        push!(V, Y_i)
        push!(I, j)
        push!(J, j)
        push!(V, Y_j)
    end
    Ybus_lcc = sparse(I, J, V, size(Ybus, 1), size(Ybus, 2))
    return Ybus_lcc
end

function _ybus_vsc(Ybus, vsc_i, vsc_j, vsc_z_i, vsc_z_j)
    I = Int[]
    J = Int[]
    V = ComplexF32[]

    for k in eachindex(vsc_i)
        i = vsc_i[k]
        j = vsc_j[k]
        z_i = vsc_z_i[k]
        z_j = vsc_z_j[k]

        Y_i = 1 / z_i
        Y_j = 1 / z_j

        push!(I, i)
        push!(J, i)
        push!(V, Y_i)
        push!(I, j)
        push!(J, j)
        push!(V, Y_j)
    end
    Ybus_vsc = sparse(I, J, V, size(Ybus, 1), size(Ybus, 2))
    return Ybus_vsc
end

function _f_lcc(
    R::Vector{Float64},
    t_i::Vector{Float64},
    t_j::Vector{Float64},
    α_i::Vector{Float64},
    α_j::Vector{Float64},
    I_dc_i::Vector{Float64},
    I_dc_j::Vector{Float64},
    x_t_i::Vector{Float64},
    x_t_j::Vector{Float64},
    lcc_P_set::Vector{Float64},
    Vm::Vector{Float64},
    i::Vector{Int64},
    j::Vector{Int64},
)
    n = length(R)
    F = Vector{Float64}(undef, 4 * n)
    for k in 1:n
        p_1 = P_CSC(t_i[k], α_i[k], I_dc_i[k], x_t_i[k], Vm[i[k]])
        p_2 = P_CSC(t_j[k], α_j[k], I_dc_j[k], x_t_j[k], Vm[j[k]])

        p_2_a = P_CSC_A(Vm[j[k]], x_t_j[k], I_dc_j[k], α_j[k]) 
        p_2_b = P_CSC_B(x_t_j[k], I_dc_j[k])
        @show t_i[k], t_j[k], α_i[k], α_j[k], I_dc_i[k], I_dc_j[k], x_t_i[k], x_t_j[k]
        @show p_1, p_2, p_2_a, p_2_b
        F1 = p_1 - lcc_P_set[k]
        F2 = p_1 + p_2 - R[k] * I_dc_i[k]^2
        F3 = α_i[k] - 0.0
        F4 = α_j[k] - π / 2
        F[4 * (k - 1) + 1] = F1
        F[4 * (k - 1) + 2] = F2
        F[4 * (k - 1) + 3] = F3
        F[4 * (k - 1) + 4] = F4
    end
    return F
end

function _legacy_J_lcc(
    pvpq::Vector{Int64},
    pq::Vector{Int64},
    Vm::Vector{Float64},
    lcc_t_i::Vector{Float64},
    lcc_t_j::Vector{Float64},
    lcc_I_dc_i::Vector{Float64},
    lcc_I_dc_j::Vector{Float64},
    lcc_α_i::Vector{Float64},
    lcc_α_j::Vector{Float64},
    lcc_x_i::Vector{Float64},
    lcc_x_j::Vector{Float64},
    size_J::Int64,
    lcc_i::Vector{Int64},
    lcc_j::Vector{Int64},
)
    # Placeholder for the actual implementation of the Jacobian for LCCs
    # This should return a SparseMatrixCSC representing the Jacobian contributions from LCCs
    I = Int[]
    J = Int[]
    V = Float64[]

    npvpq = length(pvpq)
    npq = length(pq)

    for k in 1:length(lcc_i)
        i = lcc_i[k]
        j = lcc_j[k]
        α_i = lcc_α_i[k]
        α_j = lcc_α_j[k]
        x_i = lcc_x_i[k]
        x_j = lcc_x_j[k]
        t_i = lcc_t_i[k]
        t_j = lcc_t_j[k]
        I_dc_i = lcc_I_dc_i[k]
        I_dc_j = lcc_I_dc_j[k]

        # J_C_P_D

        # J_C_Q_D

        # J_C_P_U

        i ∈ pq && push!(I, i)
        i ∈ pq && push!(J, npvpq + i)
        i ∈ pq && push!(V, t_i * sqrt(6) / π * I_dc_i * cos(α_i))

        # j ∈ pq && push!(I, j)
        # j ∈ pq && push!(J, npvpq + j)
        # j ∈ pq && push!(V, ∂P_∂u(Vm[j], t_j, -I_dc, α_j))

        # J_C_Q_U

        i ∈ pq && push!(I, npvpq + i)
        i ∈ pq && push!(J, npvpq + i)
        i ∈ pq && push!(V, ∂Q_∂u(t_i, α_i, I_dc_i, x_i, Vm[i]))

        # j ∈ pq && push!(I, npvpq + j)
        # j ∈ pq && push!(J, npvpq + j)
        # i ∈ pq && push!(V, ∂Q_∂u(t_j, α_j, -I_dc, x_j, Vm[j], θ[j]))

        # J_C_P_C

        i ∈ pvpq && push!(I, i)
        i ∈ pvpq && push!(J, npvpq + npq + 4 * (k - 1) + 1)
        i ∈ pvpq && push!(V, Vm[i] * sqrt(6) / π * I_dc_i * cos(α_i))

        i ∈ pvpq && push!(I, i)
        i ∈ pvpq && push!(J, npvpq + npq + 4 * (k - 1) + 3)
        i ∈ pvpq && push!(V, -Vm[i] * t_i * sqrt(6) / π * I_dc_i * cos(α_i) * tan(α_i))

        # j ∈ pvpq && push!(I, j)
        # j ∈ pvpq && push!(J, npvpq + npq + 4 * (k-1) + 2)
        # j ∈ pvpq && push!(V, P_CSC_A(Vm[j], t_j, -I_dc, α_j) / t_j)

        # j ∈ pvpq && push!(I, j)
        # j ∈ pvpq && push!(J, npvpq + npq + 4 * (k-1) + 4)
        # j ∈ pvpq && push!(V, -P_CSC_A(Vm[j], t_j, -I_dc, α_j) * tan(α_j))

        # J_C_Q_C

        i ∈ pq && push!(I, npvpq + i)
        i ∈ pq && push!(J, npvpq + npq + 4 * (k - 1) + 1)
        i ∈ pq && push!(V, ∂Q_∂t(t_i, α_i, I_dc_i, x_i, Vm[i]))

        i ∈ pq && push!(I, npvpq + i)
        i ∈ pq && push!(J, npvpq + npq + 4 * (k - 1) + 3)
        i ∈ pq && push!(V, ∂Q_∂α(t_i, α_i, I_dc_i, x_i, Vm[i]))

        # j ∈ pq && push!(I, npvpq + j)
        # j ∈ pq && push!(J, npvpq + npq + 4 * (k-1) + 2)
        # j ∈ pq && push!(V, ∂Q_∂t(t_j, α_j, -I_dc, x_j, Vm[j], θ[j]))

        # j ∈ pq && push!(I, npvpq + j)
        # j ∈ pq && push!(J, npvpq + npq + 4 * (k-1) + 4)
        # j ∈ pq && push!(V, ∂Q_∂α(t_j, α_j, -I_dc, x_j, Vm[j], θ[j]))

        # J_C_C_D

        # J_C_C_U

        i ∈ pq && push!(I, npvpq + npq + 4 * (k - 1) + 1)
        i ∈ pq && push!(J, npvpq + i)
        i ∈ pq && push!(V, t_i * sqrt(6) / π * I_dc_i * cos(α_i))

        i ∈ pq && push!(I, npvpq + npq + 4 * (k - 1) + 2)
        i ∈ pq && push!(J, npvpq + i)
        i ∈ pq && push!(V, t_i * sqrt(6) / π * I_dc_i * cos(α_i))

        j ∈ pq && push!(I, npvpq + npq + 4 * (k - 1) + 2)
        j ∈ pq && push!(J, npvpq + j)
        j ∈ pq && push!(V, t_j * sqrt(6) / π * I_dc_j * cos(α_j))

        # J_C_C_C

        push!(I, npvpq + npq + 4 * (k - 1) + 1)
        push!(J, npvpq + npq + 4 * (k - 1) + 1)
        push!(V, Vm[i] * sqrt(6) / π * I_dc_i * cos(α_i))

        push!(I, npvpq + npq + 4 * (k - 1) + 1)
        push!(J, npvpq + npq + 4 * (k - 1) + 3)
        push!(V, -Vm[i] * t_i * sqrt(6) / π * I_dc_i * cos(α_i) * tan(α_i))

        push!(I, npvpq + npq + 4 * (k - 1) + 2)
        push!(J, npvpq + npq + 4 * (k - 1) + 1)
        push!(V, Vm[i] * sqrt(6) / π * I_dc_i * cos(α_i))

        push!(I, npvpq + npq + 4 * (k - 1) + 2)
        push!(J, npvpq + npq + 4 * (k - 1) + 2)
        push!(V, Vm[j] * sqrt(6) / π * I_dc_j * cos(α_j))

        push!(I, npvpq + npq + 4 * (k - 1) + 2)
        push!(J, npvpq + npq + 4 * (k - 1) + 3)
        push!(V, -Vm[i] * t_i * sqrt(6) / π * I_dc_i * cos(α_i) * tan(α_i))

        push!(I, npvpq + npq + 4 * (k - 1) + 2)
        push!(J, npvpq + npq + 4 * (k - 1) + 4)
        push!(V, -Vm[j] * t_j * sqrt(6) / π * I_dc_j * cos(α_j) * tan(α_j))

        push!(I, npvpq + npq + 4 * (k - 1) + 3)
        push!(J, npvpq + npq + 4 * (k - 1) + 3)
        push!(V, 1.0)

        push!(I, npvpq + npq + 4 * (k - 1) + 4)
        push!(J, npvpq + npq + 4 * (k - 1) + 4)
        push!(V, 1.0)
    end

    return sparse(I, J, V, size_J, size_J)
end

#### LCC
###################################

ϕ(α, I_dc, x_t, t, Vm) = acos(cos(α) * sign(I_dc) - (x_t * I_dc / (sqrt(2) * Vm)))

I_ac(t, α, I_dc, x_t, Vm, Va) =
    t / Vm * sqrt(6) / π * I_dc * exp(-1im * ϕ(α, I_dc, x_t, t, Vm)) * Vm * exp(1im * Va) # check this!

Y_val(t, α, I_dc, x_t, Vm) =
    t / Vm * sqrt(6) / π * I_dc * exp(-1im * ϕ(α, I_dc, x_t, t, Vm))

Q_CSC_A(Vm, t, I_dc) = Vm * t * sqrt(6) / π * I_dc

Q_CSC_D(t, α, I_dc, x_t, Vm) = cos(α) * sign(I_dc) - x_t * I_dc / (sqrt(2) * t * Vm)

Q_CSC_C(t, α, I_dc, x_t, Vm) = acos(Q_CSC_D(t, α, I_dc, x_t, Vm))  # acos(clamp.(Q_CSC_D(t, α, I_dc, x_t, Vm), -1.0, 1.0))

Q_CSC_B(t, α, I_dc, x_t, Vm) = sin(Q_CSC_C(t, α, I_dc, x_t, Vm))

Q_CSC(t, α, I_dc, x_t, Vm) = Q_CSC_A(Vm, t, I_dc) + Q_CSC_B(t, α, I_dc, x_t, Vm)

P_CSC_A(Vm, t, I_dc, α) = Vm * t * sqrt(6) / π * I_dc * cos(α)

P_CSC_B(x_t, I_dc) = -sqrt(3 / 2) * sqrt(6) / π * x_t * I_dc^2

P_CSC(t, α, I_dc, x_t, Vm) = P_CSC_A(Vm, t, I_dc, α) + P_CSC_B(x_t, I_dc)

∂P_∂u(Vm, t, I_dc, α) = P_CSC_A(Vm, t, I_dc, α) / Vm

∂Q_∂u(t, α, I_dc, x_t, Vm) =
    Q_CSC_A(Vm, t, I_dc) * (
        Q_CSC_B(t, α, I_dc, x_t, Vm) -
        cos(Q_CSC_C(t, α, I_dc, x_t, Vm)) * x_t * I_dc /
        sqrt(1 - Q_CSC_D(t, α, I_dc, x_t, Vm)^2) * sqrt(2) * t * Vm
    ) / Vm

∂P_∂t(Vm, t, I_dc, α) = P_CSC_A(Vm, t, I_dc, α) / t

∂P_∂α(Vm, t, I_dc, α) = -P_CSC_A(Vm, t, I_dc, α) * tan(α)

∂Q_∂t(t, α, I_dc, x_t, Vm) =
    Q_CSC_A(Vm, t, I_dc) * (
        Q_CSC_B(t, α, I_dc, x_t, Vm) / t -
        cos(Q_CSC_C(t, α, I_dc, x_t, Vm)) * x_t * I_dc /
        sqrt(1 - Q_CSC_D(t, α, I_dc, x_t, Vm)^2) * sqrt(2) * t^2 * Vm
    )

∂Q_∂α(t, α, I_dc, x_t, Vm) =
    Q_CSC_A(Vm, t, I_dc) * cos(Q_CSC_C(t, α, I_dc, x_t, Vm)) *
    sin(α) * sign(I_dc) / sqrt(1 - Q_CSC_D(t, α, I_dc, x_t, Vm)^2)

∂F_∂u = P_CSC_A

∂F_∂t(Vm, t, I_dc, α) = P_CSC_A(Vm, t, I_dc, α) / t

∂F_∂α(Vm, t, I_dc, α) = -P_CSC_A(Vm, t, I_dc, α) * tan(α)

# F_CSC_1(t, α, I_dc, x_t, Vm, P_set) = P_CSC(t, α, I_dc, x_t, Vm) - P_set

# F_CSC_2(tᵢ, αᵢ, I_dc_i, I_dc_j, x_tᵢ, Vmᵢ, tⱼ, αⱼ, x_tⱼ, Vmⱼ, R) =
#     P_CSC(tᵢ, αᵢ, I_dc_i, x_tᵢ, Vmᵢ) +
#     P_CSC(tⱼ, αⱼ, I_dc_j, x_tⱼ, Vmⱼ) - R * I_dc_i^2

# μ(α, I_dc, x_t, t, Vm) = acos(cos(α) - sqrt(2) * I_dc * x_t * t / Vm) - α

# ϕ(α, μ) = acos(0.5 * cos(α) + 0.5 * cos(α + μ))

# ϕ(α, I_dc, x_t, t, Vm) = ϕ(α, μ(α, I_dc, x_t, t, Vm))

# VSC
#####################

P_VSC_ii(Vm_i, y_vsc_ii, ϕ_vsc_ii) = Vm_i^2 * y_vsc_ii * cos(ϕ_vsc_ii)
P_VSC_iq(Vm_i, Vm_q, y_vsc_iq, δ_iq, ϕ_vsc_iq) =
    Vm_i * Vm_q * y_vsc_iq * cos(ϕ_vsc_iq + δ_iq)
P_VSC(Vm_i, Vm_q, y_vsc_ii, y_vsc_iq, ϕ_vsc_ii, δ_iq, ϕ_vsc_iq) =
    P_VSC_ii(Vm_i, y_vsc_ii, ϕ_vsc_ii) + P_VSC_iq(Vm_i, Vm_q, y_vsc_iq, δ_iq, ϕ_vsc_iq)
Q_VSC_ii(Vm_i, y_vsc_ii, ϕ_vsc_ii) = Vm_i^2 * y_vsc_ii * sin(ϕ_vsc_ii)
Q_VSC_iq(Vm_i, Vm_q, y_vsc_iq, δ_iq, ϕ_vsc_iq) =
    Vm_i * Vm_q * y_vsc_iq * sin(ϕ_vsc_iq + δ_iq)
Q_VSC(Vm_i, Vm_q, y_vsc_ii, y_vsc_iq, ϕ_vsc_ii, δ_iq, ϕ_vsc_iq) =
    Q_VSC_ii(Vm_i, y_vsc_ii, ϕ_vsc_ii) + Q_VSC_iq(Vm_i, Vm_q, y_vsc_iq, δ_iq, ϕ_vsc_iq)
