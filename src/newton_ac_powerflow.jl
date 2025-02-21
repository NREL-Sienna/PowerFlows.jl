"""
    solve_powerflow!(pf::ACPowerFlow{<:ACPowerFlowSolverType}, system::PSY.System; kwargs...)

Solves the power flow in the system and writes the solution into the relevant structs.
Updates active and reactive power setpoints for generators and active and reactive
power flows for branches (calculated in the From - To direction and in the To - From direction).

Supports passing kwargs to the PF solver.

The bus types can be changed from PV to PQ if the reactive power limits are violated.

# Arguments
- `pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The power flow solver instance, can be `KLUACPowerFlow`, `HybridACPowerFlow`, or `PowerFlows.LUACPowerFlow` (to be used for testing only).
- `system::PSY.System`: The power system model.
- `kwargs...`: Additional keyword arguments.

## Keyword Arguments
- `check_connectivity::Bool`: Checks if the grid is connected. Default is `true`.
- 'check_reactive_power_limits': if `true`, the reactive power limits are enforced by changing the respective bus types from PV to PQ. Default is `false`.
- `xtol`: (only for `NLSolve`) Norm difference in `x` between two successive iterates under which convergence is declared. Default is `0.0`.
- `ftol`: (only for `NLSolve`) Infinite norm of residuals under which convergence is declared. Default is `1e-8`.
- `iterations`: (only for `NLSolve`) Maximum number of iterations. Default is `1_000`.
- `store_trace`: (only for `NLSolve`) Should a trace of the optimization algorithm's state be stored? Default is `false`.
- `show_trace`: (only for `NLSolve`) Should a trace of the optimization algorithm's state be shown on `STDOUT`? Default is `false`.
- `extended_trace`: (only for `NLSolve`) Should additional algorithm internals be added to the state trace? Default is `false`.

# Returns
- `converged::Bool`: Indicates whether the power flow solution converged.
- The power flow results are written into the system struct.

# Examples

```julia
solve_ac_powerflow!(pf, sys)

# Passing kwargs
solve_ac_powerflow!(pf, sys; check_connectivity=false)

# Passing NLsolve arguments
solve_ac_powerflow!(pf, sys; method=:newton)
```
"""
function solve_powerflow!(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    # converged must be defined in the outer scope to be visible for return
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            pf,
            system;
            check_connectivity = get(kwargs, :check_connectivity, true),
        )

        converged, V, Sbus_result =
            _ac_powerflow(data, pf, time_step; kwargs...)
        x = _calc_x(data, V, Sbus_result, time_step)

        if converged
            write_powerflow_solution!(
                system,
                x,
                data,
                get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER),
            )
            @info("PowerFlow solve converged, the results have been stored in the system")
        else
            @error("The powerflow solver returned convergence = $(converged)")
        end
    end

    return converged
end

"""
Similar to solve_powerflow!(pf, sys) but does not update the system struct with results.
Returns the results in a dictionary of dataframes.

## Examples

```julia
res = solve_powerflow(pf, sys)

# Passing NLsolve arguments
res = solve_powerflow(pf, sys; method=:newton)
```
"""
function solve_powerflow(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    system::PSY.System;
    kwargs...,
)
    # df_results must be defined in the oueter scope first to be visible for return
    df_results = Dict{String, DataFrames.DataFrame}()
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(
            pf,
            system;
            check_connectivity = get(kwargs, :check_connectivity, true),
        )

        converged, V, Sbus_result = _ac_powerflow(data, pf, time_step; kwargs...)
        x = _calc_x(data, V, Sbus_result, time_step)

        if converged
            @info("PowerFlow solve converged, the results are exported in DataFrames")
            df_results = write_results(pf, system, data, x, time_step)
        else
            df_results = missing
            @error("The powerflow solver returned convergence = $(converged)")
        end
    end

    return df_results
end

"""
    solve_powerflow!(data::ACPowerFlowData; pf::ACPowerFlow{<:ACPowerFlowSolverType} = ACPowerFlow(), kwargs...)

Solve the multiperiod AC power flow problem for the given power flow data.

The bus types can be changed from PV to PQ if the reactive power limits are violated.

# Arguments
- `data::ACPowerFlowData`: The power flow data containing netwthe grid information and initial conditions.
- `pf::ACPowerFlow{<:ACPowerFlowSolverType}`: The power flow solver type. Defaults to `KLUACPowerFlow`.
- `kwargs...`: Additional keyword arguments.

# Keyword Arguments
- `check_connectivity::Bool`: Checks if the grid is connected. Default is `true`.
- 'check_reactive_power_limits': if `true`, the reactive power limits are enforced by changing the respective bus types from PV to PQ. Default is `false`.
- `time_steps`: Specifies the time steps to solve. Defaults to sorting and collecting the keys of `data.timestep_map`.

# Returns
- `Nothing`. The results are written directly to the `data` object.

# Description
This function solves the AC power flow problem for each time step specified in `data`. 
It preallocates memory for the results and iterates over the sorted time steps. 
    For each time step, it calls the `_ac_powerflow` function to solve the power flow equations and updates the `data` object with the results. 
    If the power flow converges, it updates the active and reactive power injections, as well as the voltage magnitudes and angles for different bus types (REF, PV, PQ). 
    If the power flow does not converge, it sets the corresponding entries in `data` to `NaN`. 
    Finally, it calculates the branch power flows and updates the `data` object.

# Notes
- If the grid topology changes (e.g., tap positions of transformers or in-service status of branches), the admittance matrices `Yft` and `Ytf` must be updated.
- If `Yft` and `Ytf` change between time steps, the branch flow calculations must be moved inside the loop.

# Examples
```julia
solve_powerflow!(data)
```
"""
function solve_powerflow!(
    data::ACPowerFlowData;
    pf::ACPowerFlow{<:ACPowerFlowSolverType} = ACPowerFlow(),
    enable_progress_bar::Bool = true,
    kwargs...,
)
    sorted_time_steps = get(kwargs, :time_steps, sort(collect(keys(data.timestep_map))))
    # preallocate results
    ts_converged = fill(false, length(sorted_time_steps))
    ts_V = zeros(Complex{Float64}, first(size(data.bus_type)), length(sorted_time_steps))
    ts_S = zeros(Complex{Float64}, first(size(data.bus_type)), length(sorted_time_steps))

    # TODO If anything in the grid topology changes, 
    #  e.g. tap positions of transformers or in service 
    #  status of branches, Yft and Ytf must be updated!
    Yft = data.power_network_matrix.yft
    Ytf = data.power_network_matrix.ytf
    fb = data.power_network_matrix.fb
    tb = data.power_network_matrix.tb

    progress_bar = ProgressMeter.Progress(
        length(sorted_time_steps);
        enabled = enable_progress_bar,
        desc = "Multi-period power flow",
        showspeed = true,
    )

    for time_step in sorted_time_steps
        ProgressMeter.update!(
            progress_bar,
            time_step;
            showvalues = [
                (:Step, time_step),
            ],
        )
        converged, V, Sbus_result =
            _ac_powerflow(data, pf, time_step; kwargs...)
        ts_converged[time_step] = converged
        ts_V[:, time_step] .= V
        ts_S[:, time_step] .= Sbus_result

        if converged
            ref, pv, pq = bus_type_idx(data, time_step)

            # write results to data:
            # write results for REF
            data.bus_activepower_injection[ref, time_step] .=
                real.(Sbus_result[ref]) .+ data.bus_activepower_withdrawals[ref, time_step]
            data.bus_reactivepower_injection[ref, time_step] .=
                imag.(Sbus_result[ref]) .+
                data.bus_reactivepower_withdrawals[ref, time_step]
            # write Q results for PV
            data.bus_reactivepower_injection[pv, time_step] .=
                imag.(Sbus_result[pv]) .+ data.bus_reactivepower_withdrawals[pv, time_step]
            # results for PQ buses do not need to be updated -> already consistent with inputs

            # write voltage results
            data.bus_magnitude[pq, time_step] .= abs.(V[pq])
            data.bus_angles[pq, time_step] .= angle.(V[pq])
            data.bus_angles[pv, time_step] .= angle.(V[pv])
        else
            data.bus_activepower_injection[:, time_step] .= NaN
            data.bus_activepower_withdrawals[:, time_step] .= NaN
            data.bus_reactivepower_injection[:, time_step] .= NaN
            data.bus_reactivepower_withdrawals[:, time_step] .= NaN
            data.bus_magnitude[:, time_step] .= NaN
            data.bus_angles[:, time_step] .= NaN
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

function _ac_powerflow(
    data::ACPowerFlowData,
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    time_step::Int64;
    kwargs...,
)
    check_reactive_power_limits = get(kwargs, :check_reactive_power_limits, false)

    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged, V, Sbus_result =
            _newton_powerflow(pf, data, time_step; kwargs...)
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
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
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
    check_reactive_power_limits,
    time_step::Int64;
    kwargs...,
)
    for _ in 1:MAX_REACTIVE_POWER_ITERATIONS
        converged, V, Sbus_result =
            _newton_powerflow(pf, data, time_step; kwargs...)
        if !converged || !check_reactive_power_limits ||
           _check_q_limit_bounds!(data, Sbus_result, time_step)
            return converged, V, Sbus_result
        end
    end
    # todo: throw error? set converged to false? -> leave as is for now
    @error("could not enforce reactive power limits after $MAX_REACTIVE_POWER_ITERATIONS")
    return converged, V, Sbus_result
end

function bus_type_idx(
    data::ACPowerFlowData,
    time_step::Int64 = 1,
    bus_types::Tuple{Vararg{PSY.ACBusTypes}} = (
        PSY.ACBusTypes.REF,
        PSY.ACBusTypes.PV,
        PSY.ACBusTypes.PQ,
    ),
)
    # Find indices for each bus type
    return [
        findall(x -> x == bus_type, data.bus_type[:, time_step]) for bus_type in bus_types
    ]
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

function _update_dSbus_dV!(rows::Vector{Int32}, cols::Vector{Int32},
    V::Vector{Complex{Float64}}, Ybus::SparseMatrixCSC{Complex{Float64}, Int64},
    diagV::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    diagVnorm::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    diagIbus::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    diagIbus_diag::Vector{Complex{Float64}},
    dSbus_dVa::SparseMatrixCSC{Complex{Float64}, Int32},
    dSbus_dVm::SparseMatrixCSC{Complex{Float64}, Int32},
    r_dSbus_dVa::SparseMatrixCSC{Float64, Int32},
    r_dSbus_dVm::SparseMatrixCSC{Float64, Int32},
    i_dSbus_dVa::SparseMatrixCSC{Float64, Int32},
    i_dSbus_dVm::SparseMatrixCSC{Float64, Int32},
    Ybus_diagVnorm::SparseMatrixCSC{Complex{Float64}, Int32},
    conj_Ybus_diagVnorm::SparseMatrixCSC{Complex{Float64}, Int32},
    diagV_conj_Ybus_diagVnorm::SparseMatrixCSC{Complex{Float64}, Int32},
    conj_diagIbus::LinearAlgebra.Diagonal{Complex{Float64}, Vector{Complex{Float64}}},
    conj_diagIbus_diagVnorm::LinearAlgebra.Diagonal{
        Complex{Float64},
        Vector{Complex{Float64}},
    },
    Ybus_diagV::SparseMatrixCSC{Complex{Float64}, Int32},
    conj_Ybus_diagV::SparseMatrixCSC{Complex{Float64}, Int32})
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
    Sbus_result::Vector{Complex{Float64}},
    time_step::Int64,
)
    Vm = abs.(V)
    Va = angle.(V)
    n_buses = length(V)
    # mock the expected x format, where the values depend on the type of the bus:
    x = zeros(Float64, 2 * n_buses)
    bus_types = view(data.bus_type, :, time_step)
    for (ix, bt) in enumerate(bus_types)
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
    rows::Vector{Int32},
    cols::Vector{Int32},
    pvpq::Vector{Int64},
    pq::Vector{Int64},
)
    J_block = SparseMatrixCSC{Float64, Int32}(
        sparse(rows, cols, 0.0, maximum(rows), maximum(cols)),
    )
    J = [J_block[pvpq, pvpq] J_block[pvpq, pq]; J_block[pq, pvpq] J_block[pq, pq]]
    return J
end

function _newton_powerflow(
    pf::ACPowerFlow{KLUACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...,
)
    # Fetch maxIterations and tol from kwargs, or use defaults if not provided
    maxIterations = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    tol = get(kwargs, :tol, DEFAULT_NR_TOL)
    i = 0

    Ybus = data.power_network_matrix.data

    # Find indices for each bus type
    ref, pv, pq = bus_type_idx(data, time_step)
    pvpq = [pv; pq]

    # nref = length(ref)
    npv = length(pv)
    npq = length(pq)
    npvpq = npv + npq

    Vm = data.bus_magnitude[:, time_step]
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

    # Convert rows and cols to Int32
    rows = Int32.(rows)
    cols = Int32.(cols)

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
    cache = KLULinSolveCache(J)

    while i < maxIterations && !converged
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
            symbolic_factor!(cache, J)
        end

        try
            # factorize the numeric object of KLU inplace, while reusing the symbolic object
            numeric_refactor!(cache, J)

            # solve inplace - the results are written to F, so that we must use F instead of dx for updating V
            solve!(cache, F)
        catch e
            @error("KLU factorization failed: $e")
            return (converged, V, Sbus_result)
        end

        # KLU.solve! overwrites F with the solution instead of returning it as dx, so -F is used here to update V
        _update_V!(F, V, Vm, Va, pv, pq, dx_Va_pv, dx_Va_pq, dx_Vm_pq)

        # here F is mismatch again
        _update_F!(F, Sbus_result, mis, dx_Va_pv, dx_Va_pq, dx_Vm_pq, V, Ybus, Sbus, pv, pq)

        converged = norm(F, Inf) < tol
    end

    if !converged
        V .*= NaN
        Sbus_result .*= NaN
        @error("The powerflow solver with KLU did not converge after $i iterations")
    else
        if pf.calc_loss_factors
            data.loss_factors[ref, :] .= 0.0
            penalty_factors!(
                J,
                collect(real.(hcat(dSbus_dVa[ref, pvpq], dSbus_dVm[ref, pq]))[:]),
                view(data.loss_factors, pvpq, time_step), collect(1:npvpq))
        end
        @info("The powerflow solver with KLU converged after $i iterations")
    end

    return (converged, V, Sbus_result)
end
