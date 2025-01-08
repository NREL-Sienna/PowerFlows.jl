"""
Solves a the power flow into the system and writes the solution into the relevant structs.
Updates generators active and reactive power setpoints and branches active and reactive
power flows (calculated in the From - To direction) (see
[`flow_val`](@ref))

Supports passing NLsolve kwargs in the args. By default shows the solver trace.

Arguments available for `nlsolve`:

  - `get_connectivity::Bool`: Checks if the network is connected. Default true
  - `method` : See NLSolve.jl documentation for available solvers
  - `xtol`: norm difference in `x` between two successive iterates under which
    convergence is declared. Default: `0.0`.
  - `ftol`: infinite norm of residuals under which convergence is declared.
    Default: `1e-8`.
  - `iterations`: maximum number of iterations. Default: `1_000`.
  - `store_trace`: should a trace of the optimization algorithm's state be
    stored? Default: `false`.
  - `show_trace`: should a trace of the optimization algorithm's state be shown
    on `STDOUT`? Default: `false`.
  - `extended_trace`: should additifonal algorithm internals be added to the state
    trace? Default: `false`.

## Examples

```julia
solve_ac_powerflow!(sys)
# Passing NLsolve arguments
solve_ac_powerflow!(sys, method=:newton)
```
"""
function solve_ac_powerflow!(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    check_reactive_power_limits = get(kwargs, :check_reactive_power_limits, false)
    data = PowerFlowData(
        pf,
        system;
        check_connectivity = get(kwargs, :check_connectivity, true),
    )
    max_iterations = DEFAULT_MAX_REDISTRIBUTION_ITERATIONS
    converged, x = _solve_powerflow!(pf, data, check_reactive_power_limits; kwargs...)
    if converged
        write_powerflow_solution!(system, x, max_iterations)
        @info("PowerFlow solve converged, the results have been stored in the system")
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return converged
    end
    @error("The powerflow solver returned convergence = $(converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return converged
end

"""
Similar to solve_powerflow!(sys) but does not update the system struct with results.
Returns the results in a dictionary of dataframes.

## Examples

```julia
res = solve_powerflow(sys)
# Passing NLsolve arguments
res = solve_powerflow(sys, method=:newton)
```
"""
function solve_powerflow(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    data = PowerFlowData(
        pf,
        system;
        check_connectivity = get(kwargs, :check_connectivity, true),
    )

    converged, x = _solve_powerflow!(pf, data, pf.check_reactive_power_limits; kwargs...)

    if converged
        @info("PowerFlow solve converged, the results are exported in DataFrames")
        df_results = write_results(pf, system, data, x)
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return df_results
    end
    @error("The powerflow solver returned convergence = $(converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return converged
end

function _check_q_limit_bounds!(data::ACPowerFlowData, zero::Vector{Float64})
    bus_names = data.power_network_matrix.axes[1]
    within_limits = true
    for (ix, b) in enumerate(data.bus_type)
        if b == PSY.ACBusTypes.PV
            Q_gen = zero[2 * ix - 1]
        else
            continue
        end

        if Q_gen <= data.bus_reactivepower_bounds[ix][1]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix] = data.bus_reactivepower_bounds[ix][1]
        elseif Q_gen >= data.bus_reactivepower_bounds[ix][2]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix] = data.bus_reactivepower_bounds[ix][2]
        else
            @debug "Within Limits"
        end
    end
    return within_limits
end

function _solve_powerflow!(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    data::ACPowerFlowData,
    check_reactive_power_limits;
    nlsolve_kwargs...,
)
    if check_reactive_power_limits
        for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
            converged, x = _newton_powerflow(pf, data; nlsolve_kwargs...)
            if converged
                if _check_q_limit_bounds!(data, x)
                    return converged, x
                end
            else
                return converged, x
            end
        end
    else
        return _newton_powerflow(pf, data; nlsolve_kwargs...)
    end
end

function _newton_powerflow(
    pf::ACPowerFlow{NLSolveACPowerFlow},
    data::ACPowerFlowData;
    nlsolve_kwargs...,
)
    pf = PolarPowerFlow(data)
    J = PowerFlows.PolarPowerFlowJacobian(data, pf.x0)

    df = NLsolve.OnceDifferentiable(pf, J, pf.x0, pf.residual, J.Jv)
    res = NLsolve.nlsolve(df, pf.x0; nlsolve_kwargs...)
    if !res.f_converged
        @error(
            "The powerflow solver NLSolve did not converge (returned convergence = $(res.f_converged))"
        )
    end
    return res.f_converged, res.zero
end

function _update_V!(dx::Vector{Float64}, V::Vector{Complex{Float64}}, Vm::Vector{Float64},
    Va::Vector{Float64},
    pv::Vector{Int64}, pq::Vector{Int64}, dx_Va_pv::Vector{Int64}, dx_Va_pq::Vector{Int64},
    dx_Vm_pq::Vector{Int64})
    for (i, j) in zip(pv, dx_Va_pv)
        Va[i] -= dx[j]
    end

    for (i, j) in zip(pq, dx_Va_pq)
        Va[i] -= dx[j]
    end

    for (i, j) in zip(pq, dx_Vm_pq)
        Vm[i] -= dx[j]
    end

    V .= Vm .* exp.(1im .* Va)

    Vm .= abs.(V)
    Va .= angle.(V)
    return
end

function _update_F!(F::Vector{Float64}, mis::Vector{Complex{Float64}},
    dx_Va_pv::Vector{Int64}, dx_Va_pq::Vector{Int64}, dx_Vm_pq::Vector{Int64},
    V::Vector{Complex{Float64}}, Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
    Sbus::Vector{Complex{Float64}},
    pv::Vector{Int64}, pq::Vector{Int64})

    #mis .= V .* conj.(Ybus * V) .- Sbus
    LinearAlgebra.mul!(mis, Ybus, V)
    mis .= V .* conj.(mis) .- Sbus

    for (i, j) in zip(dx_Va_pv, pv)
        F[i] = real(mis[j])
    end

    for (i, j) in zip(dx_Va_pq, pq)
        F[i] = real(mis[j])
    end

    for (i, j) in zip(dx_Vm_pq, pq)
        F[i] = imag(mis[j])
    end
    return
end

function _update_dSbus_dV!(rows::Vector{Int64}, cols::Vector{Int64},
    V::Vector{Complex{Float64}}, Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
    diagV::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    diagVnorm::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    diagIbus::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    diagIbus_diag::Vector{Complex{Float64}},
    dSbus_dVa::SparseMatrixCSC{Complex{Float64}, Int64},
    dSbus_dVm::SparseMatrixCSC{Complex{Float64}, Int64},
    r_dSbus_dVa::SparseMatrixCSC{Float64, Int64},
    r_dSbus_dVm::SparseMatrixCSC{Float64, Int64},
    i_dSbus_dVa::SparseMatrixCSC{Float64, Int64},
    i_dSbus_dVm::SparseMatrixCSC{Float64, Int64},
    Ybus_diagVnorm::SparseMatrixCSC{Complex{Float64}, Int64},
    conj_Ybus_diagVnorm::SparseMatrixCSC{Complex{Float64}, Int64},
    diagV_conj_Ybus_diagVnorm::SparseMatrixCSC{Complex{Float64}, Int64},
    conj_diagIbus::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    conj_diagIbus_diagVnorm::LinearAlgebra.Diagonal{
        Complex{Float64},
        Vector{Complex{Float64}},
    },
    Ybus_diagV::SparseMatrixCSC{Complex{Float64}, Int64},
    conj_Ybus_diagV::SparseMatrixCSC{Complex{Float64}, Int64})
    for i in eachindex(V)
        diagV[i, i] = V[i]
        diagVnorm[i, i] = V[i] / abs(V[i])
    end

    # manually calculate the diagIbus matrix
    LinearAlgebra.mul!(diagIbus_diag, Ybus, V)
    for i in eachindex(V)
        diagIbus[i, i] = diagIbus_diag[i]
    end

    # use the available matrices temporarily to calculate the dSbus_dV matrices
    # original formula:
    # dSbus_dVm .= diagV * conj.(Ybus * diagVnorm) + conj.(diagIbus) * diagVnorm
    # non-allocating version:

    LinearAlgebra.mul!(Ybus_diagVnorm, Ybus, diagVnorm)
    conj_Ybus_diagVnorm .= conj.(Ybus_diagVnorm)
    LinearAlgebra.mul!(diagV_conj_Ybus_diagVnorm, diagV, conj_Ybus_diagVnorm)
    conj_diagIbus .= conj.(diagIbus)
    LinearAlgebra.mul!(conj_diagIbus_diagVnorm, conj_diagIbus, diagVnorm)

    dSbus_dVm .= diagV_conj_Ybus_diagVnorm
    LinearAlgebra.axpy!(1, conj_diagIbus_diagVnorm, dSbus_dVm)

    # original formula:
    # dSbus_dVa .= 1im * diagV * conj.(diagIbus - Ybus * diagV)
    # non-allocating version:
    LinearAlgebra.mul!(Ybus_diagV, Ybus, diagV)
    # Take the conjugate of the result (conj(Ybus * diagV)); conj_diagIbus is already available
    conj_Ybus_diagV .= conj.(Ybus_diagV)

    # write the result of conj.(diagIbus - Ybus * diagV) in conj_Ybus_diagV
    LinearAlgebra.axpby!(1, conj_diagIbus, -1, conj_Ybus_diagV)

    # Multiply the result by diagV
    LinearAlgebra.mul!(dSbus_dVa, diagV, conj_Ybus_diagV)

    # Now multiply by 1im to get the final result
    #dSbus_dVa .*= 1im
    LinearAlgebra.mul!(dSbus_dVa, dSbus_dVa, 1im)

    # this loop is slower so we should use vectorize assignments below
    # for c in cols
    #     for r in rows
    #         r_dSbus_dVa[r, c] = real(dSbus_dVa[r, c])
    #         i_dSbus_dVa[r, c] = imag(dSbus_dVa[r, c])
    #         r_dSbus_dVm[r, c] = real(dSbus_dVm[r, c])
    #         i_dSbus_dVm[r, c] = imag(dSbus_dVm[r, c])
    #     end
    # end

    # sometimes can allocate so we have to use the for loop above
    r_dSbus_dVa .= real.(dSbus_dVa)
    r_dSbus_dVm .= real.(dSbus_dVm)
    i_dSbus_dVa .= imag.(dSbus_dVa)
    i_dSbus_dVm .= imag.(dSbus_dVm)
    return
end

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

function _update_submatrix!(
    A::SparseMatrixCSC,
    B::SparseMatrixCSC,
    rows_A::Vector{Int64},
    cols_A::Vector{Int64},
    rows_B::Vector{Int64},
    cols_B::Vector{Int64},
)
    for idj in eachindex(cols_A)
        for idi in eachindex(rows_A)
            A[rows_A[idi], cols_A[idj]] = B[rows_B[idi], cols_B[idj]]
        end
    end
    return
end

function _update_J!(J::SparseMatrixCSC,
    r_dSbus_dVa::SparseMatrixCSC,
    r_dSbus_dVm::SparseMatrixCSC,
    i_dSbus_dVa::SparseMatrixCSC,
    i_dSbus_dVm::SparseMatrixCSC,
    pvpq::Vector{Int64},
    pq::Vector{Int64},
    j_pvpq::Vector{Int64},
    j_pq::Vector{Int64},
)
    _update_submatrix!(J, r_dSbus_dVa, j_pvpq, j_pvpq, pvpq, pvpq)
    _update_submatrix!(J, r_dSbus_dVm, j_pvpq, j_pq, pvpq, pq)
    _update_submatrix!(J, i_dSbus_dVa, j_pq, j_pvpq, pq, pvpq)
    _update_submatrix!(J, i_dSbus_dVm, j_pq, j_pq, pq, pq)
    return
end

function _calc_x(
    data::ACPowerFlowData,
    V::Vector{Complex{Float64}},
    Va::Vector{Float64},
    Vm::Vector{Float64},
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
    n_buses::Int64,
)
    # mock the expected x format, where the values depend on the type of the bus:
    x = zeros(Float64, 2 * n_buses)
    Sbus_result = V .* conj(Ybus * V)  # todo preallocate and pass as parameter
    for (ix, b) in enumerate(data.bus_type)
        if b == PSY.ACBusTypes.REF
            # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
            x[2 * ix - 1] = real(Sbus_result[ix]) + data.bus_activepower_withdrawals[ix]
            x[2 * ix] = imag(Sbus_result[ix]) + data.bus_reactivepower_withdrawals[ix]
        elseif b == PSY.ACBusTypes.PV
            # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
            x[2 * ix - 1] = imag(Sbus_result[ix]) + data.bus_reactivepower_withdrawals[ix]
            x[2 * ix] = Va[ix]
        elseif b == PSY.ACBusTypes.PQ
            # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
            x[2 * ix - 1] = Vm[ix]
            x[2 * ix] = Va[ix]
        end
    end
    return x
end

function _newton_powerflow(
    pf::ACPowerFlow{KLUACPowerFlow},
    data::ACPowerFlowData;
    nlsolve_kwargs...,
)
    # Fetch maxIter and tol from kwargs, or use defaults if not provided
    maxIter = get(nlsolve_kwargs, :maxIter, DEFAULT_NR_MAX_ITER)
    tol = get(nlsolve_kwargs, :tol, DEFAULT_NR_TOL)
    i = 0

    Ybus = data.power_network_matrix.data

    # Find indices for each bus type
    ref = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.REF, data.bus_type)
    pv = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PV, data.bus_type)
    pq = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PQ, data.bus_type)
    pvpq = [pv; pq]

    # nref = length(ref)
    npv = length(pv)
    npq = length(pq)
    npvpq = npv + npq
    n_buses = length(data.bus_type)

    Vm = data.bus_magnitude[:]
    # prevent unfeasible starting values for Vm; for pv and ref buses we cannot do this:
    Vm[pq] .= clamp.(Vm[pq], 0.9, 1.1)
    Va = data.bus_angles[:]
    V = zeros(Complex{Float64}, length(Vm))
    V .= Vm .* exp.(1im .* Va)

    # early return if only ref buses present - no need to solve the power flow
    converged = npvpq == 0
    if converged
        # if only ref buses present, we do not need to enter the power flow loop
        x = _calc_x(data, V, Va, Vm, Ybus, n_buses)
        return (converged, x)
    end

    # we need to define lookups for mappings of pv, pq buses onto the internal J indexing
    pvpq_lookup = zeros(Int64, maximum([ref; pvpq]) + 1)
    pvpq_lookup[pvpq] .= 1:npvpq
    pq_lookup = zeros(Int64, maximum([ref; pvpq]) + 1)
    pq_lookup[pq] .= 1:npq

    # define the internal J indexing using the lookup arrays
    j_pvpq = pvpq_lookup[pvpq]
    j_pq = npvpq .+ pq_lookup[pq]

    # indices for updating of V
    dx_Va_pv = Vector{Int64}([1:npv...])
    dx_Va_pq = Vector{Int64}([(npv + 1):(npv + npq)...])
    dx_Vm_pq = Vector{Int64}([(npv + npq + 1):(npv + 2 * npq)...])

    Sbus =
        data.bus_activepower_injection[:] - data.bus_activepower_withdrawals[:] +
        1im * (data.bus_reactivepower_injection[:] - data.bus_reactivepower_withdrawals[:])

    # Pre-allocate mis and F and create views for the respective real and imaginary sections of the arrays:
    mis = zeros(Complex{Float64}, length(V))
    F = zeros(Float64, npvpq + npq)
    mis .= V .* conj.(Ybus * V) .- Sbus
    F .= [real(mis[pvpq]); imag(mis[pq])]

    # preallocate Jacobian matrix and arrays for calculating dSbus_dVa, dSbus_dVm
    rows, cols = SparseArrays.findnz(Ybus)

    #diagV = sparse(1:n_buses, 1:n_buses, V)
    diagV = LinearAlgebra.Diagonal(V)
    diagIbus_diag = zeros(Complex{Float64}, size(V, 1))
    diagIbus = LinearAlgebra.Diagonal(diagIbus_diag)

    diagVnorm = LinearAlgebra.Diagonal(V ./ abs.(V))

    Ybus_diagVnorm = sparse(rows, cols, Complex{Float64}(0))
    conj_Ybus_diagVnorm = sparse(rows, cols, Complex{Float64}(0))
    diagV_conj_Ybus_diagVnorm = sparse(rows, cols, Complex{Float64}(0))
    conj_diagIbus = conj.(diagIbus)
    conj_diagIbus_diagVnorm = conj.(diagIbus)
    Ybus_diagV = sparse(rows, cols, Complex{Float64}(0))
    conj_Ybus_diagV = sparse(rows, cols, Complex{Float64}(0))

    dSbus_dVm = sparse(rows, cols, Complex{Float64}(0))
    dSbus_dVa = sparse(rows, cols, Complex{Float64}(0))
    r_dSbus_dVa = sparse(rows, cols, Float64(0))
    r_dSbus_dVm = sparse(rows, cols, Float64(0))
    i_dSbus_dVa = sparse(rows, cols, Float64(0))
    i_dSbus_dVm = sparse(rows, cols, Float64(0))

    # maybe use this in the future?
    # pvpq_rows = pvpq_lookup[rows][pvpq_lookup[rows] .!= 0]
    # pvpq_cols = pvpq_lookup[cols][pvpq_lookup[cols] .!= 0]
    # pq_rows = pq_lookup[rows][pq_lookup[rows] .!= 0] 
    # pq_cols = pq_lookup[cols][pq_lookup[cols] .!= 0]

    J_block = sparse(rows, cols, Float64(0), maximum(rows), maximum(cols))
    J = [J_block[pvpq, pvpq] J_block[pvpq, pq]; J_block[pq, pvpq] J_block[pq, pq]]

    # preallocate the KLU factorization object - symbolic object only
    colptr = KLU.decrement(J.colptr)
    rowval = KLU.decrement(J.rowval)
    n = size(J, 1)
    factor_J = KLU.KLUFactorization(n, colptr, rowval, J.nzval)
    KLU.klu_analyze!(factor_J)
    rf = Ref(factor_J.common)
    # factorization for the numeric object does not work here:
    # factor_J._numeric = KLU.klu_l_factor(colptr, rowval, J.nzval, factor_J._symbolic, rf)

    while i < maxIter && !converged
        i += 1

        _update_dSbus_dV!(rows, cols, V, Ybus, diagV, diagVnorm, diagIbus, diagIbus_diag,
            dSbus_dVa, dSbus_dVm, r_dSbus_dVa, r_dSbus_dVm, i_dSbus_dVa, i_dSbus_dVm,
            Ybus_diagVnorm, conj_Ybus_diagVnorm, diagV_conj_Ybus_diagVnorm,
            conj_diagIbus, conj_diagIbus_diagVnorm, Ybus_diagV, conj_Ybus_diagV)

        # todo: improve pvpq, pq, j_pvpq, j_pq (use more specific indices)
        _update_J!(
            J,
            r_dSbus_dVa,
            r_dSbus_dVm,
            i_dSbus_dVa,
            i_dSbus_dVm,
            pvpq,
            pq,
            j_pvpq,
            j_pq,
        )

        # Workaround for the issue with KLU.klu_l_factor
        # background: KLU.klu_l_factor does not work properly with the preallocated J matrix with dummy values
        # the workaround is to initialize the numeric object here in the loop once and then refactorize the matrix in the loop inplace
        if i == 1
            # works when J values are ok:
            factor_J._numeric =
                KLU.klu_l_factor(colptr, rowval, J.nzval, factor_J._symbolic, rf)
        end

        # factorize the numeric object of KLU inplace, while reusing the symbolic object
        KLU.klu_l_refactor(
            colptr,
            rowval,
            J.nzval,
            factor_J._symbolic,
            factor_J._numeric,
            rf,
        )

        # solve inplace - the results are written to F, so that we must use F instead of dx for updating V
        KLU.klu_l_solve(
            factor_J._symbolic,
            factor_J._numeric,
            size(F, 1),
            size(F, 2),
            F,
            rf,
        )

        # KLU.solve! overwrites F with the solution instead of returning it as dx, so -F is used here to update V
        _update_V!(F, V, Vm, Va, pv, pq, dx_Va_pv, dx_Va_pq, dx_Vm_pq)

        # here F is mismatch again
        _update_F!(F, mis, dx_Va_pv, dx_Va_pq, dx_Vm_pq, V, Ybus, Sbus, pv, pq)

        converged = LinearAlgebra.norm(F, Inf) < tol
    end

    if !converged
        @error("The powerflow solver with KLU did not converge after $i iterations")
    else
        @info("The powerflow solver with KLU converged after $i iterations")
    end

    x = _calc_x(data, V, Va, Vm, Ybus, n_buses)

    return (converged, x)
end

# legacy NR implementation - here we do not care about allocations, we use this function only for testing purposes
function _newton_powerflow(
    pf::ACPowerFlow{LUACPowerFlow},
    data::ACPowerFlowData;
    nlsolve_kwargs...,
)
    # Fetch maxIter and tol from kwargs, or use defaults if not provided
    maxIter = get(nlsolve_kwargs, :maxIter, DEFAULT_NR_MAX_ITER)
    tol = get(nlsolve_kwargs, :tol, DEFAULT_NR_TOL)
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

    Vm = data.bus_magnitude[:]
    # prevent unfeasible starting values for Vm; for pv and ref buses we cannot do this:
    Vm[pq] .= clamp.(Vm[pq], 0.9, 1.1)
    Va = data.bus_angles[:]
    V = zeros(Complex{Float64}, length(Vm))
    V .= Vm .* exp.(1im .* Va)

    # pre-allocate dx
    dx = zeros(Float64, npv + 2 * npq)

    Ybus = data.power_network_matrix.data

    Sbus =
        data.bus_activepower_injection[:] - data.bus_activepower_withdrawals[:] +
        1im * (data.bus_reactivepower_injection[:] - data.bus_reactivepower_withdrawals[:])

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

    x = _calc_x(data, V, Va, Vm, Ybus, n_buses)

    return (converged, x)
end
