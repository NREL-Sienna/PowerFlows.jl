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

function _newton_powerflow(
    pf::ACPowerFlow{KLUACPowerFlow},
    data::ACPowerFlowData;
    nlsolve_kwargs...,
)
    # Fetch maxIter and tol from kwargs, or use defaults if not provided
    maxIter = get(nlsolve_kwargs, :maxIter, DEFAULT_NR_MAX_ITER)
    tol = get(nlsolve_kwargs, :tol, DEFAULT_NR_TOL)
    i = 0

    # Find indices for each bus type
    ref = findall(x -> x == PowerSystems.ACBusTypesModule.ACBusTypes.REF, data.bus_type)
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
    @. Vm[pq] = clamp.(Vm[pq], 0.9, 1.1)
    Va = data.bus_angles[:]
    V = zeros(Complex{Float64}, length(Vm))
    @. V = Vm .* exp.(1im * Va)

    Va_pv = view(Va, pv)
    Va_pq = view(Va, pq)
    Vm_pq = view(Vm, pq)

    # pre-allocate dx
    dx = zeros(Float64, npv + 2 * npq)

    dx_Va_pv = view(dx, 1:npv)
    dx_Va_pq = view(dx, (npv + 1):(npv + npq))
    dx_Vm_pq = view(dx, (npv + npq + 1):(npv + 2 * npq))

    Ybus = data.power_network_matrix.data

    Sbus =
        data.bus_activepower_injection[:] - data.bus_activepower_withdrawals[:] +
        1im * (data.bus_reactivepower_injection[:] - data.bus_reactivepower_withdrawals[:])
    
    # Pre-allocate mis and F and create views for the respective real and imaginary sections of the arrays:
    mis = zeros(Complex{Float64}, length(V))
    mis_pvpq = view(mis, pvpq)
    mis_pq = view(mis, pq)

    F = zeros(Float64, npvpq + npq)
    F_real = view(F, 1:npvpq)
    F_imag = view(F, npvpq + 1:npvpq + npq)

    mis .= V .* conj.(Ybus * V) .- Sbus
    @. F_real = real(mis_pvpq)  # In-place assignment to the real part, using views
    @. F_imag = imag(mis_pq)  # In-place assignment to the imaginary part, using views

    converged = npvpq == 0  # if only ref buses present, we do not need to enter the loop

    # preallocate Jacobian matrix
    rows = vcat(1:npvpq, 1:npvpq, npvpq+1:npvpq+npq, npvpq+1:npvpq+npq)
    cols = vcat(1:npvpq, npvpq+1:npvpq+npq, 1:npvpq, npvpq+1:npvpq+npq)
    J = sparse(rows, cols, Float64(0))

    # we need to define lookups for mappings of pv, pq buses onto the internal J indexing
    pvpq_lookup = zeros(Int64, maximum([ref; pvpq]) + 1)
    pvpq_lookup[pvpq] .= 1:npvpq
    pq_lookup = zeros(Int64, maximum([ref; pvpq]) + 1)
    pq_lookup[pq] .= 1:npq

    # with the pre-allocated J and lookups, we can define views into the sub-matrices of the J matrix for updating the J matrix in the NR loop
    j11 = view(J, pvpq_lookup[pvpq], pvpq_lookup[pvpq])
    j12 = view(J, pvpq_lookup[pvpq], npvpq .+ pq_lookup[pq])
    j21 = view(J, npvpq .+ pq_lookup[pq], pvpq_lookup[pvpq])
    j22 = view(J, npvpq .+ pq_lookup[pq], npvpq .+ pq_lookup[pq])
    
    # we need views of the diagonals to avoid using LinearAlgebra.Diagonal:
    diagV = sparse(1:n_buses, 1:n_buses, Complex{Float64}(1))
    diag_idx = LinearAlgebra.diagind(diagV)
    diagV_diag = view(diagV, diag_idx)

    diagIbus = sparse(1:n_buses, 1:n_buses, Complex{Float64}(1))
    diagIbus_diag = view(diagIbus, diag_idx)

    diagVnorm = sparse(1:n_buses, 1:n_buses, Complex{Float64}(1))
    diagVnorm_diag = view(diagVnorm, diag_idx)

    # pre-allocate the dSbus_dVm, dSbus_dVa to have the same structure as Ybus 
    # they will follow the structure of Ybus except maybe when Ybus has zero values in its diagonal, which we do not expect here
    #rows, cols, _ = SparseArrays.findnz(Ybus)
    #dSbus_dVm = sparse(rows, cols, Complex{Float64}(0))
    #dSbus_dVa = sparse(rows, cols, Complex{Float64}(0))

    # preallocate dSbus_dVm, dSbus_dVa with correct structure:
    dSbus_dVm = diagV * conj.(Ybus * diagVnorm) + conj.(diagIbus) * diagVnorm
    dSbus_dVa = 1im * diagV * conj.(diagIbus - Ybus * diagV)

    # create views for the sub-arrays of Sbus_dVa, Sbus_dVm for updating the J:
    Sbus_dVa_j11 = view(dSbus_dVa, pvpq, pvpq)
    Sbus_dVm_j12 = view(dSbus_dVm, pvpq, pq)
    Sbus_dVa_j21 = view(dSbus_dVa, pq, pvpq)
    Sbus_dVm_j22 = view(dSbus_dVm, pq, pq)

    while i < maxIter && !converged
        i += 1
 
        ## use the new value of V to update dSbus_dVa, dSbus_dVm:
        diagV_diag .= V
        diagIbus_diag .= Ybus * V
        @. diagVnorm_diag = V ./ abs.(V)
        dSbus_dVm .= diagV * conj.(Ybus * diagVnorm) + conj.(diagIbus) * diagVnorm
        dSbus_dVa .= 1im * diagV * conj.(diagIbus - Ybus * diagV)

        # update the Jacobian by setting values through the pre-defined views for j11, j12, j21, j22
        @. j11 = real(Sbus_dVa_j11)
        @. j12 = real(Sbus_dVm_j12)
        @. j21 = imag(Sbus_dVa_j21)
        @. j22 = imag(Sbus_dVm_j22)

        factor_J = KLU.klu(J)
        dx .= -(factor_J \ F)

        Va_pv .+= dx_Va_pv
        Va_pq .+= dx_Va_pq
        Vm_pq .+= dx_Vm_pq

        @. V = Vm .* exp.(1im * Va)

        Vm .= abs.(V)
        Va .= angle.(V)

        mis .= V .* conj.(Ybus * V) .- Sbus
        @. F_real = real(mis_pvpq)  # In-place assignment to the real part
        @. F_imag = imag(mis_pq)  # In-place assignment to the imaginary part
        converged = LinearAlgebra.norm(F, Inf) < tol
    end

    if !converged
        @error("The powerflow solver with KLU did not converge after $i iterations")
    else
        @info("The powerflow solver with KLU converged after $i iterations")
    end

    # mock the expected x format, where the values depend on the type of the bus:
    n_buses = length(data.bus_type)
    x = zeros(Float64, 2 * n_buses)
    Sbus_result = V .* conj(Ybus * V)
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

    return converged, x
end
