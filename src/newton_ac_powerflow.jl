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
function solve_powerflow!(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    time_step::Int64 = 1,
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

    converged, V, Sbus_result = _ac_powereflow(data, pf; time_step = time_step, kwargs...)
    x = _calc_x(data, V, Sbus_result)

    if converged
        write_powerflow_solution!(
            system,
            x,
            data,
            get(kwargs, :maxIter, DEFAULT_NR_MAX_ITER),
        )
        @info("PowerFlow solve converged, the results have been stored in the system")
    else
        @error("The powerflow solver returned convergence = $(converged)")
    end

    #Restore original per unit base
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

    converged, V, Sbus_result = _ac_powereflow(data, pf; kwargs...)
    x = _calc_x(data, V, Sbus_result)

    if converged
        @info("PowerFlow solve converged, the results are exported in DataFrames")
        df_results = write_results(pf, system, data, x)
    else
        df_results = missing
        @error("The powerflow solver returned convergence = $(converged)")
    end

    #Restore original per unit base
    PSY.set_units_base_system!(system, settings_unit_cache)

    return df_results
end

# Multiperiod power flow - work in progress
function solve_powerflow(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    data::ACPowerFlowData,
    system::PSY.System;
    kwargs...,
)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")

    sorted_time_steps = sort(collect(keys(data.timestep_map)))
    # preallocate results
    ts_converged = zeros(Bool, 1, length(sorted_time_steps))
    ts_V = zeros(Complex{Float64}, length(data.bus_type[:, 1]), length(sorted_time_steps))
    ts_S = zeros(Complex{Float64}, length(data.bus_type[:, 1]), length(sorted_time_steps))

    results = Dict()

    for t in sorted_time_steps
        converged, V, Sbus_result =
            _ac_powereflow(data, pf; time_step = t, kwargs...)
        ts_converged[1, t] = converged
        ts_V[:, t] .= V
        ts_S[:, t] .= Sbus_result

        # temporary implementation that will need to be improved:
        if converged
            x = _calc_x(data, V, Sbus_result)
            results[data.timestep_map[t]] = write_results(pf, system, data, x)
        else
            results[data.timestep_map[t]] = missing
        end

        #todo: implement write_results for multiperiod power flow

        # if converged
        #     @info("PowerFlow solve converged, the results are exported in DataFrames")
        #     df_results = write_results(pf, system, data, x)
        # else
        #     df_results = missing
        #     @error("The powerflow solver returned convergence = $(converged)")
        # end
    end

    #Restore original per unit base
    PSY.set_units_base_system!(system, settings_unit_cache)

    # todo:
    # return df_results

    return results
end

# Multiperiod power flow - work in progress
function solve_powerflow!(
    data::ACPowerFlowData;
    kwargs...,
)
    pf = ACPowerFlow()  # todo: somehow store in data which PF to use (see issue #50)

    sorted_time_steps = get(kwargs, :time_steps, sort(collect(keys(data.timestep_map))))
    # preallocate results
    ts_converged = fill(false, length(sorted_time_steps))
    ts_V = zeros(Complex{Float64}, length(data.bus_type[:, 1]), length(sorted_time_steps))
    ts_S = zeros(Complex{Float64}, length(data.bus_type[:, 1]), length(sorted_time_steps))

    # TODO If anything in the grid topology changes, 
    #  e.g. tap positions of transformers or in service 
    #  status of branches, Yft and Ytf must be updated!
    Yft = data.power_network_matrix.yft
    Ytf = data.power_network_matrix.ytf
    fb = data.power_network_matrix.fb
    tb = data.power_network_matrix.tb

    for t in sorted_time_steps
        converged, V, Sbus_result =
            _ac_powereflow(data, pf; time_step = t, kwargs...)
        ts_converged[t] = converged
        ts_V[:, t] .= V
        ts_S[:, t] .= Sbus_result

        if converged
            ref = findall(
                x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.REF,
                data.bus_type[:, t],
            )
            pv = findall(
                x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PV,
                data.bus_type[:, t],
            )
            pq = findall(
                x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PQ,
                data.bus_type[:, t],
            )

            # temporary implementation that will need to be improved:
            # write results for REF
            data.bus_activepower_injection[ref, t] .=
                real.(Sbus_result[ref]) .+ data.bus_activepower_withdrawals[ref, t]
            data.bus_reactivepower_injection[ref, t] .=
                imag.(Sbus_result[ref]) .+ data.bus_reactivepower_withdrawals[ref, t]
            # write Q results for PV
            data.bus_reactivepower_injection[pv, t] .=
                imag.(Sbus_result[pv]) .+ data.bus_reactivepower_withdrawals[pv, t]
            # results for PQ buses do not need to be updated -> already consistent with inputs

            # write voltage results
            data.bus_magnitude[pq, t] .= abs.(V[pq])
            data.bus_angles[pq, t] .= angle.(V[pq])
            data.bus_angles[pv, t] .= angle.(V[pv])
        else
            data.bus_activepower_injection[:, t] .= NaN64
            data.bus_activepower_withdrawals[:, t] .= NaN64
            data.bus_reactivepower_injection[:, t] .= NaN64
            data.bus_reactivepower_withdrawals[:, t] .= NaN64
            data.bus_magnitude[:, t] .= NaN64
            data.bus_angles[:, t] .= NaN64
        end
    end

    # write branch flows
    # TODO if Yft, Ytf change between time steps, this must be moved inside the loop!
    Sft = ts_V[fb, :] .* conj.(Yft * ts_V)
    Stf = ts_V[tb, :] .* conj.(Ytf * ts_V)

    data.branch_activepower_flow_from_to .= real.(Sft)
    data.branch_reactivepower_flow_from_to .= imag.(Sft)
    data.branch_activepower_flow_to_from .= real.(Stf)
    data.branch_reactivepower_flow_to_from .= imag.(Stf)

    data.converged .= ts_converged

    return
end

function _ac_powereflow(
    data::ACPowerFlowData,
    pf::ACPowerFlow{<:ACPowerFlowSolverType};
    time_step::Int64 = 1,
    kwargs...,
)
    check_reactive_power_limits = get(kwargs, :check_reactive_power_limits, false)

    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged, V, Sbus_result =
            _newton_powerflow(pf, data; time_step = time_step, kwargs...)
        if !converged || !check_reactive_power_limits ||
           _check_q_limit_bounds!(data, Sbus_result, time_step)
            return (converged, V, Sbus_result)
        end
    end

    @error("could not enforce reactive power limits after $MAX_REACTIVE_POWER_ITERATIONS")
    return (converged, V, Sbus_result)
end

function _check_q_limit_bounds!(
    data::ACPowerFlowData,
    Sbus_result::Vector{Complex{Float64}},
    time_step::Int64,
)
    bus_names = data.power_network_matrix.axes[1]
    within_limits = true
    for (ix, bt) in enumerate(data.bus_type[:, time_step])
        if bt == PSY.ACBusTypes.PV
            Q_gen = imag(Sbus_result[ix])
        else
            continue
        end

        if Q_gen <= data.bus_reactivepower_bounds[ix, time_step][1]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix, time_step] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix, time_step] =
                data.bus_reactivepower_bounds[ix, time_step][1]
        elseif Q_gen >= data.bus_reactivepower_bounds[ix, time_step][2]
            @info "Bus $(bus_names[ix]) changed to PSY.ACBusTypes.PQ"
            within_limits = false
            data.bus_type[ix, time_step] = PSY.ACBusTypes.PQ
            data.bus_reactivepower_injection[ix, time_step] =
                data.bus_reactivepower_bounds[ix, time_step][2]
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
    time_step::Int64 = 1,
    kwargs...,
)
    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged, V, Sbus_result =
            _newton_powerflow(pf, data; time_step = time_step, kwargs...)
        if !converged || !check_reactive_power_limits ||
           _check_q_limit_bounds!(data, Sbus_result, time_step)
            return converged, V, Sbus_result
        end
    end
    # todo: throw error? set converged to false?
    @error("could not enforce reactive power limits after $MAX_REACTIVE_POWER_ITERATIONS")
    return converged, V, Sbus_result
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

function _update_F!(F::Vector{Float64}, Sbus_result::Vector{Complex{Float64}},
    mis::Vector{Complex{Float64}},
    dx_Va_pv::Vector{Int64}, dx_Va_pq::Vector{Int64}, dx_Vm_pq::Vector{Int64},
    V::Vector{Complex{Float64}}, Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
    Sbus::Vector{Complex{Float64}},
    pv::Vector{Int64}, pq::Vector{Int64})

    #mis .= V .* conj.(Ybus * V) .- Sbus
    LinearAlgebra.mul!(Sbus_result, Ybus, V)
    Sbus_result .= V .* conj(Sbus_result)
    mis .= Sbus_result .- Sbus

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
    Sbus_result::Vector{Complex{Float64}};
    time_step::Int64 = 1,
)
    Vm = abs.(V)
    Va = angle.(V)
    n_buses = length(V)
    # mock the expected x format, where the values depend on the type of the bus:
    x = zeros(Float64, 2 * n_buses)
    for (ix, bt) in enumerate(data.bus_type[:, time_step])
        if bt == PSY.ACBusTypes.REF
            # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
            x[2 * ix - 1] = real(Sbus_result[ix]) + data.bus_activepower_withdrawals[ix]
            x[2 * ix] = imag(Sbus_result[ix]) + data.bus_reactivepower_withdrawals[ix]
        elseif bt == PSY.ACBusTypes.PV
            # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
            x[2 * ix - 1] = imag(Sbus_result[ix]) + data.bus_reactivepower_withdrawals[ix]
            x[2 * ix] = Va[ix]
        elseif bt == PSY.ACBusTypes.PQ
            # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
            x[2 * ix - 1] = Vm[ix]
            x[2 * ix] = Va[ix]
        end
    end
    return x
end

function _preallocate_J(
    rows::Vector{Int64},
    cols::Vector{Int64},
    pvpq::Vector{Int64},
    pq::Vector{Int64},
)
    J_block = sparse(rows, cols, Float64(0), maximum(rows), maximum(cols))
    J = [J_block[pvpq, pvpq] J_block[pvpq, pq]; J_block[pq, pvpq] J_block[pq, pq]]
    return J
end

function _newton_powerflow(
    pf::ACPowerFlow{KLUACPowerFlow},
    data::ACPowerFlowData;
    time_step::Int64 = 1,
    kwargs...,
)
    # Fetch maxIter and tol from kwargs, or use defaults if not provided
    maxIter = get(kwargs, :maxIter, DEFAULT_NR_MAX_ITER)
    tol = get(kwargs, :tol, DEFAULT_NR_TOL)
    i = 0

    solver_data = data.solver_data[time_step]

    Ybus = data.power_network_matrix.data

    # Find indices for each bus type
    ref = findall(
        x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.REF,
        data.bus_type[:, time_step],
    )
    pv = findall(
        x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PV,
        data.bus_type[:, time_step],
    )
    pq = findall(
        x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.PQ,
        data.bus_type[:, time_step],
    )
    pvpq = [pv; pq]

    # nref = length(ref)
    npv = length(pv)
    npq = length(pq)
    npvpq = npv + npq

    Vm = data.bus_magnitude[:, time_step]
    # prevent unfeasible starting values for Vm; for pv and ref buses we cannot do this:
    Vm[pq] .= clamp.(Vm[pq], 0.9, 1.1)
    Va = data.bus_angles[:, time_step]
    V = zeros(Complex{Float64}, length(Vm))
    V .= Vm .* exp.(1im .* Va)

    # early return if only ref buses present - no need to solve the power flow
    converged = npvpq == 0
    if converged
        # if only ref buses present, we do not need to enter the power flow loop
        Sbus_result = V .* conj(Ybus * V)
        return (converged, V, Sbus_result)
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
        data.bus_activepower_injection[:, time_step] -
        data.bus_activepower_withdrawals[:, time_step] +
        1im * (
            data.bus_reactivepower_injection[:, time_step] -
            data.bus_reactivepower_withdrawals[:, time_step]
        )

    # Pre-allocate mis and F and create views for the respective real and imaginary sections of the arrays:
    mis = zeros(Complex{Float64}, length(V))
    Sbus_result = zeros(Complex{Float64}, length(V))
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

    J = _preallocate_J(rows, cols, pvpq, pq)

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

        try
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

        catch e
            @error("KLU factorization failed: $e")
            return (converged, V, Sbus_result)
        end

        # KLU.solve! overwrites F with the solution instead of returning it as dx, so -F is used here to update V
        _update_V!(F, V, Vm, Va, pv, pq, dx_Va_pv, dx_Va_pq, dx_Vm_pq)

        # here F is mismatch again
        _update_F!(F, Sbus_result, mis, dx_Va_pv, dx_Va_pq, dx_Vm_pq, V, Ybus, Sbus, pv, pq)

        converged = LinearAlgebra.norm(F, Inf) < tol
    end

    if !converged
        V .*= NaN64
        Sbus_result .*= NaN64
        @error("The powerflow solver with KLU did not converge after $i iterations")
    else
        solver_data.J = J
        solver_data.dSbus_dV_ref =
            [vec(real.(dSbus_dVa[ref, :][:, pvpq])); vec(real.(dSbus_dVm[ref, :][:, pq]))]
        @info("The powerflow solver with KLU converged after $i iterations")
    end

    return (converged, V, Sbus_result)
end
