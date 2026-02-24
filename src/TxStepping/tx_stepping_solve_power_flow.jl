"""Default lambda schedule for the TxStepping homotopy method.
At lambda=1 the system is trivially solvable; at lambda=0 the original system is recovered."""
const TX_STEPPING_LAMBDA_SCHEDULE = [1.0, 0.9, 0.7, 0.5, 0.3, 0.1, 0.0]

"""
    _ref_mag_for_bus(data, time_step) -> Vector{Float64}

Return a vector of length n where entry `i` is the voltage magnitude of the
REF bus in bus `i`'s subnetwork. Used to interpolate PV voltage setpoints
toward V_ref at lambda=1 and to initialize the flat start.
"""
function _ref_mag_for_bus(data::PowerFlowData, time_step::Int)
    n = size(data.bus_type, 1)
    bus_lookup = data.power_network_matrix.lookup[1]
    # subnetwork_axes is keyed by REF bus number (set by assign_reference_buses!).
    subnetwork_axes = data.power_network_matrix.subnetwork_axes

    ref_mag = zeros(Float64, n)
    for (ref_bus_num, (sub_buses, _)) in subnetwork_axes
        ref_vm = data.bus_magnitude[bus_lookup[ref_bus_num], time_step]
        for bus_num in sub_buses
            ref_mag[bus_lookup[bus_num]] = ref_vm
        end
    end
    return ref_mag
end
const TX_STEPPING_GAMMA = 1e4
const TX_STEPPING_DEFAULT_MAX_ITER = 100
const TX_STEPPING_DEFAULT_TOL = 1e-9

"""
Run Newton-Raphson inner loop for a single lambda step of TxStepping.
Returns `true` if converged.
"""
function _tx_stepping_newton!(
    state::TxSteppingSolverState,
    residual::TxSteppingResidual,
    jacobian::TxSteppingJacobian,
    linSolveCache::KLULinSolveCache{J_INDEX_TYPE},
    time_step::Int,
    lambda::Float64;
    maxIterations::Int = TX_STEPPING_DEFAULT_MAX_ITER,
    tol::Float64 = TX_STEPPING_DEFAULT_TOL,
)
    data = residual.data

    # compute initial residual
    compute_residual!(residual, time_step, lambda)
    if norm(residual.Rv, Inf) < tol
        return true
    end

    n = size(data.bus_type, 1)
    for iter in 1:maxIterations
        # update Jacobian
        _update_jacobian!(
            jacobian.Jv, jacobian.rhs, state, data,
            jacobian.ref_mag, time_step, lambda,
        )

        # factor and solve: Jv * x_new = rhs
        # _update_jacobian! sets up the system so that the solution is the new state,
        # not a delta (ref bus rows enforce V_new = V_set directly).
        numeric_refactor!(linSolveCache, jacobian.Jv)
        solve!(linSolveCache, jacobian.rhs)
        x_new = jacobian.rhs  # solution is written in-place

        # update state with backtracking line search
        new_V = reinterpret(ComplexF64, @view(x_new[1:(2 * n)]))
        dV = new_V .- state.V
        n_pv = length(x_new) - 2 * n
        dq = zeros(Float64, n)
        if n_pv > 0
            pv_ix = 1
            bus_types_view = @view data.bus_type[:, time_step]
            for (ix, bt) in enumerate(bus_types_view)
                bt != PSY.ACBusTypes.PV && continue
                dq[ix] = x_new[2 * n + pv_ix] - state.q_g[ix]
                pv_ix += 1
            end
        end

        r_prev = norm(residual.Rv, Inf)
        V_prev = copy(state.V)
        q_g_prev = copy(state.q_g)
        alpha = 1.0
        for _ls in 1:10
            state.V .= V_prev .+ alpha .* dV
            state.q_g .= q_g_prev .+ alpha .* dq
            compute_residual!(residual, time_step, lambda)
            if norm(residual.Rv, Inf) < r_prev || alpha <= 1.0 / 512
                break
            end
            alpha *= 0.5
        end

        if norm(residual.Rv, Inf) < tol
            @info "TxStepping NR converged in $iter iterations"
            return true
        end
    end
    @warn "TxStepping NR failed to converge after $maxIterations iterations " *
          "(residual = $(norm(residual.Rv, Inf)))"
    return false
end

"""
    solve_power_flow!(data::TxSteppingPowerFlowData; kwargs...)

Solve the TxStepping homotopy power flow for the given data.

# Keyword Arguments
- `tol::Float64`: convergence tolerance. Default: $TX_STEPPING_DEFAULT_TOL
- `maxIterations::Int`: max NR iterations per lambda step. Default: $TX_STEPPING_DEFAULT_MAX_ITER
- `gamma::Float64`: homotopy scaling parameter. Default: $TX_STEPPING_GAMMA
- `lambda_schedule::Vector{Float64}`: homotopy steps from 1 to 0. Default: $TX_STEPPING_LAMBDA_SCHEDULE
"""
function solve_power_flow!(
    data::TxSteppingPowerFlowData;
    kwargs...,
)
    pf = get_pf(data)
    merged_kwargs = merge(get_solver_kwargs(pf), kwargs)
    tol = get(merged_kwargs, :tol, TX_STEPPING_DEFAULT_TOL)
    maxIterations = get(merged_kwargs, :maxIterations, TX_STEPPING_DEFAULT_MAX_ITER)
    gamma = get(merged_kwargs, :gamma, TX_STEPPING_GAMMA)
    lambda_schedule = get(merged_kwargs, :lambda_schedule, TX_STEPPING_LAMBDA_SCHEDULE)
    time_step = 1  # single time step for now

    n = size(data.bus_type, 1)
    bus_types = @view data.bus_type[:, time_step]

    # Precompute per-bus REF magnitude for Vset interpolation.
    ref_mag = _ref_mag_for_bus(data, time_step)

    # Initialize state: at lambda=1, all buses start at their subnetwork's V_ref.
    V0 = ref_mag .* ones(ComplexF64, n)
    # Set REF bus voltages exactly (they include the angle).
    for (ix, bt) in enumerate(bus_types)
        bt != PSY.ACBusTypes.REF && continue
        V0[ix] =
            data.bus_magnitude[ix, time_step] *
            exp(im * data.bus_angles[ix, time_step])
    end
    q_g0 = zeros(Float64, n)
    state = TxSteppingSolverState(V0, q_g0)

    n_pv = count(==(PSY.ACBusTypes.PV), bus_types)
    residual = TxSteppingResidual(data, state, zeros(Float64, 2 * n + n_pv), ref_mag)
    jacobian = _build_tx_stepping_jacobian(data, state, ref_mag, time_step)

    # symbolic factorization (done once, reused across all lambda steps)
    linSolveCache = KLULinSolveCache(jacobian.Jv)
    symbolic_factor!(linSolveCache, jacobian.Jv)

    converged = false
    for lambda in lambda_schedule
        @info "TxStepping: lambda = $lambda"
        PNM.update_y_lambda!(data.power_network_matrix, lambda, gamma)

        # At intermediate lambda steps the residual floor is set by gamma * eps,
        # so use a relaxed tolerance; only enforce the user tolerance at lambda=0.
        step_tol = lambda > 0 ? max(tol, 1e-6) : tol
        converged = _tx_stepping_newton!(
            state, residual, jacobian, linSolveCache, time_step, lambda;
            maxIterations = maxIterations, tol = step_tol,
        )

        if !converged
            @error "TxStepping failed to converge at lambda = $lambda"
            data.converged[time_step] = false
            return false
        end
    end

    # Write converged voltages back to data
    for ix in 1:n
        data.bus_magnitude[ix, time_step] = abs(state.V[ix])
        data.bus_angles[ix, time_step] = angle(state.V[ix])
    end

    # Update reactive power injections at PV buses
    for (ix, bt) in enumerate(bus_types)
        bt != PSY.ACBusTypes.PV && continue
        data.bus_reactive_power_injections[ix, time_step] += state.q_g[ix]
    end

    # Update active and reactive power injections at REF buses from converged voltages.
    # S_net = V * conj(Y * V), then P_gen = P_net + P_withdrawals, Q_gen = Q_net + Q_withdrawals.
    # y_lambda has no off-diagonal entries in REF rows (for Jacobian efficiency), so
    # reconstruct I_ref = sum_j Y_ij * V_j using the transpose: Y_ij = Y_ji for
    # symmetric Ybus, but we use I_ref = Y_diag * V_ref + sum_{j!=ref} Y_ji^T * V_j.
    # Simpler: compute from y_series + y_shunt which have the full data.
    y_series = data.power_network_matrix.y_series
    y_shunt = data.power_network_matrix.y_shunt
    for (ix, bt) in enumerate(bus_types)
        bt != PSY.ACBusTypes.REF && continue
        # I_ref = (y_series[ref,:] + y_shunt[ref]*e_ref) * V
        # y_series is symmetric, so row ref = column ref.
        I_ref = y_shunt[ix] * state.V[ix]
        for k in SparseArrays.nzrange(y_series, ix)
            j = SparseArrays.rowvals(y_series)[k]
            I_ref += SparseArrays.nonzeros(y_series)[k] * state.V[j]
        end
        S_net = state.V[ix] * conj(I_ref)
        data.bus_active_power_injections[ix, time_step] =
            real(S_net) + data.bus_active_power_withdrawals[ix, time_step]
        data.bus_reactive_power_injections[ix, time_step] =
            imag(S_net) + data.bus_reactive_power_withdrawals[ix, time_step]
    end

    data.converged[time_step] = true
    return true
end

# --- Public API ---

"""
    solve_power_flow(pf::TxSteppingPowerFlow, system::PSY.System; kwargs...)

Solve the TxStepping homotopy power flow and return results as a dictionary of DataFrames.
"""
function solve_power_flow(
    pf::TxSteppingPowerFlow,
    system::PSY.System;
    kwargs...,
)
    df_results = Dict{String, DataFrames.DataFrame}()
    converged = false
    time_step = 1
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(pf, system)
        converged = solve_power_flow!(data; kwargs...)
        if converged
            @info(
                "TxStepping PowerFlow solve converged, the results are exported in DataFrames"
            )
            df_results = write_results(pf, system, data, time_step)
        else
            df_results = missing
            @error("The TxStepping power flow solver returned convergence = $(converged)")
        end
    end
    return df_results
end

"""
    solve_and_store_power_flow!(pf::TxSteppingPowerFlow, system::PSY.System; kwargs...)

Solve the TxStepping homotopy power flow and write the solution into the system structs.
"""
function solve_and_store_power_flow!(
    pf::TxSteppingPowerFlow,
    system::PSY.System;
    kwargs...,
)
    converged = false
    with_units_base(system, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(pf, system)
        converged = solve_power_flow!(data; kwargs...)
        if converged
            write_power_flow_solution!(
                system,
                pf,
                data,
                get(kwargs, :maxIterations, TX_STEPPING_DEFAULT_MAX_ITER),
            )
            @info(
                "TxStepping PowerFlow solve converged, the results have been stored in the system"
            )
        else
            @error("The TxStepping power flow solver returned convergence = $converged")
        end
    end
    return converged
end
