
"""
Solves a the power flow into the system and writes the solution into the relevant structs.
Updates generators active and reactive power setpoints and branches active and reactive
power flows (calculated in the From - To direction) (see
[`flow_val`](@ref))

Supports solving using Finite Differences Method (instead of using analytic Jacobian)
by setting finite_diff = true.
Supports passing NLsolve kwargs in the args. By default shows the solver trace.

Arguments available for `nlsolve`:

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
run_powerflow!(sys)
# Passing NLsolve arguments
run_powerflow!(sys, method=:newton)
# Using Finite Differences
run_powerflow!(sys, finite_diff=true)
```
"""
function run_powerflow!(system::PSY.System; finite_diff=false, kwargs...)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    res = _run_powerflow(system, finite_diff; kwargs...)
    if res.f_converged
        write_powerflow_solution!(system, res.zero)
        @info("PowerFlow solve converged, the results have been stored in the system")
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return res.f_converged
    end
    @error("The powerflow solver returned convergence = $(res.f_converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return res.f_converged
end

"""
Similar to run_powerflow!(sys) but does not update the system struct with results.
Returns the results in a dictionary of dataframes.

## Examples

```julia
res = run_powerflow(sys)
# Passing NLsolve arguments
res = run_powerflow(sys, method=:newton)
# Using Finite Differences
res = run_powerflow(sys, finite_diff=true)
```
"""
function run_powerflow(system::PSY.System; finite_diff=false, kwargs...)
    #Save per-unit flag
    settings_unit_cache = deepcopy(system.units_settings.unit_system)
    #Work in System per unit
    PSY.set_units_base_system!(system, "SYSTEM_BASE")
    res = _run_powerflow(system, finite_diff; kwargs...)
    if res.f_converged
        @info("PowerFlow solve converged, the results are exported in DataFrames")
        df_results = write_results(system, res.zero)
        #Restore original per unit base
        PSY.set_units_base_system!(system, settings_unit_cache)
        return df_results
    end
    @error("The powerflow solver returned convergence = $(res.f_converged)")
    PSY.set_units_base_system!(system, settings_unit_cache)
    return res.f_converged
end

function _run_powerflow(system::PSY.System, finite_diff::Bool; kwargs...)
    ybus_kw = [:check_connectivity, :connectivity_method]
    ybus_kwargs = (k for k in kwargs if first(k) ∈ ybus_kw)
    kwargs = (k for k in kwargs if first(k) ∉ ybus_kw)

    buses = sort!(collect(PSY.get_components(PSY.Bus, system)), by=x -> PSY.get_number(x))
    N_BUS = length(buses)

    # assumes the ordering in YPSY.Bus is the same as in the buses.
    Yb = PSY.Ybus(system; ybus_kwargs...).data
    a = collect(1:N_BUS)
    I, J, V = SparseArrays.findnz(Yb)
    neighbors = [Set{Int}([i]) for i in 1:N_BUS]
    for nz in eachindex(V)
        push!(neighbors[I[nz]], J[nz])
        push!(neighbors[J[nz]], I[nz])
    end
    x0 = zeros(N_BUS * 2)

    #Create Jacobian structure
    J0_I = Int[]
    J0_J = Int[]
    J0_V = Float64[]

    for ix_f in a
        F_ix_f_r = 2 * ix_f - 1
        F_ix_f_i = 2 * ix_f

        for ix_t in neighbors[ix_f]
            X_ix_t_fst = 2 * ix_t - 1
            X_ix_t_snd = 2 * ix_t
            b = PSY.get_bustype(buses[ix_t])
            #Set to 0.0 only on connected buses
            if b == PSY.BusTypes.REF
                if ix_f == ix_t
                    #Active PF w/r Local Active Power
                    push!(J0_I, F_ix_f_r)
                    push!(J0_J, X_ix_t_fst)
                    push!(J0_V, 0.0)
                    #Rective PF w/r Local Reactive Power
                    push!(J0_I, F_ix_f_i)
                    push!(J0_J, X_ix_t_snd)
                    push!(J0_V, 0.0)
                end
            elseif b == PSY.BusTypes.PV
                #Active PF w/r Angle
                push!(J0_I, F_ix_f_r)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
                #Reactive PF w/r Angle
                push!(J0_I, F_ix_f_i)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
                if ix_f == ix_t
                    #Reactive PF w/r Local Reactive Power
                    push!(J0_I, F_ix_f_i)
                    push!(J0_J, X_ix_t_fst)
                    push!(J0_V, 0.0)
                end
            elseif b == PSY.BusTypes.PQ
                #Active PF w/r VoltageMag
                push!(J0_I, F_ix_f_r)
                push!(J0_J, X_ix_t_fst)
                push!(J0_V, 0.0)
                #Reactive PF w/r VoltageMag
                push!(J0_I, F_ix_f_i)
                push!(J0_J, X_ix_t_fst)
                push!(J0_V, 0.0)
                #Active PF w/r Angle
                push!(J0_I, F_ix_f_r)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
                #Reactive PF w/r Angle
                push!(J0_I, F_ix_f_i)
                push!(J0_J, X_ix_t_snd)
                push!(J0_V, 0.0)
            end
        end
    end
    J0 = SparseArrays.sparse(J0_I, J0_J, J0_V)

    # Use vectors to cache data for closure
    # These should be read only
    P_GEN_BUS = fill(0.0, N_BUS)
    Q_GEN_BUS = fill(0.0, N_BUS)
    P_LOAD_BUS = fill(0.0, N_BUS)
    Q_LOAD_BUS = fill(0.0, N_BUS)
    P_net = fill(0.0, N_BUS)
    Q_net = fill(0.0, N_BUS)
    Vm = fill(0.0, N_BUS)
    θ = fill(0.0, N_BUS)

    state_variable_count = 1
    sources =
        PSY.get_components(PSY.StaticInjection, system, d -> !isa(d, PSY.ElectricLoad))

    for (ix, b) in enumerate(buses)
        bus_angle = PSY.get_angle(b)::Float64
        Vm[ix] = bus_voltage_magnitude = PSY.get_magnitude(b)::Float64
        P_GEN_BUS[ix] = 0.0
        Q_GEN_BUS[ix] = 0.0
        PSY.get_ext(b)["neighbors"] = neighbors[ix]
        for gen in sources
            !(PSY.get_available(gen) && PSY.get_status(gen)) && continue
            if gen.bus == b
                P_GEN_BUS[ix] += PSY.get_active_power(gen)
                Q_GEN_BUS[ix] += PSY.get_reactive_power(gen)
            end
        end

        P_LOAD_BUS[ix], Q_LOAD_BUS[ix] = _get_load_data(system, b)
        if b.bustype == PSY.BusTypes.REF
            injection_components = PSY.get_components(
                PSY.StaticInjection,
                system,
                d -> is_available_source(d, b),
            )
            isempty(injection_components) && throw(
                IS.ConflictingInputsError(
                    "The slack bus does not have any injection component. Power Flow can not proceed",
                ),
            )
            θ[ix] = PSY.get_angle(b)::Float64
            x0[state_variable_count] = P_GEN_BUS[ix]
            x0[state_variable_count + 1] = Q_GEN_BUS[ix]
            state_variable_count += 2
        elseif b.bustype == PSY.BusTypes.PV
            x0[state_variable_count] = Q_GEN_BUS[ix]
            x0[state_variable_count + 1] = bus_angle
            state_variable_count += 2
        elseif b.bustype == PSY.BusTypes.PQ
            x0[state_variable_count] = bus_voltage_magnitude
            x0[state_variable_count + 1] = bus_angle
            state_variable_count += 2
        else
            throw(ArgumentError("PSY.Bustype not recognized"))
        end
    end

    @assert state_variable_count - 1 == N_BUS * 2
    bus_types = PSY.get_bustype.(buses)
    function pf!(F::Vector{Float64}, X::Vector{Float64})
        for (ix, b) in enumerate(bus_types)
            if b == PSY.BusTypes.REF
                # When bustype == REFERENCE PSY.Bus, state variables are Active and Reactive Power Generated
                P_net[ix] = X[2 * ix - 1] - P_LOAD_BUS[ix]
                Q_net[ix] = X[2 * ix] - Q_LOAD_BUS[ix]
            elseif b == PSY.BusTypes.PV
                # When bustype == PV PSY.Bus, state variables are Reactive Power Generated and Voltage Angle
                P_net[ix] = P_GEN_BUS[ix] - P_LOAD_BUS[ix]
                Q_net[ix] = X[2 * ix - 1] - Q_LOAD_BUS[ix]
                θ[ix] = X[2 * ix]
            elseif b == PSY.BusTypes.PQ
                # When bustype == PQ PSY.Bus, state variables are Voltage Magnitude and Voltage Angle
                P_net[ix] = P_GEN_BUS[ix] - P_LOAD_BUS[ix]
                Q_net[ix] = Q_GEN_BUS[ix] - Q_LOAD_BUS[ix]
                Vm[ix] = X[2 * ix - 1]
                θ[ix] = X[2 * ix]
            end
        end

        # F is active and reactive power balance equations at all buses
        for ix_f in a
            S_re = -P_net[ix_f]
            S_im = -Q_net[ix_f]
            for ix_t in neighbors[ix_f]
                gb = real(Yb[ix_f, ix_t])
                bb = imag(Yb[ix_f, ix_t])
                if ix_f == ix_t
                    S_re += Vm[ix_f] * Vm[ix_t] * gb
                    S_im += -Vm[ix_f] * Vm[ix_t] * bb
                else
                    S_re +=
                        Vm[ix_f] *
                        Vm[ix_t] *
                        (gb * cos(θ[ix_f] - θ[ix_t]) + bb * sin(θ[ix_f] - θ[ix_t]))
                    S_im +=
                        Vm[ix_f] *
                        Vm[ix_t] *
                        (gb * sin(θ[ix_f] - θ[ix_t]) - bb * cos(θ[ix_f] - θ[ix_t]))
                end
            end
            F[2 * ix_f - 1] = S_re
            F[2 * ix_f] = S_im
        end
    end

    function jsp!(J::SparseArrays.SparseMatrixCSC{Float64, Int}, X::Vector{Float64})
        for ix_f in a
            F_ix_f_r = 2 * ix_f - 1
            F_ix_f_i = 2 * ix_f

            for ix_t in neighbors[ix_f]
                X_ix_t_fst = 2 * ix_t - 1
                X_ix_t_snd = 2 * ix_t
                b = bus_types[ix_t]

                if b == PSY.BusTypes.REF
                    # State variables are Active and Reactive Power Generated
                    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - x[2*i-1]
                    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
                    # x does not appears in p_flow and q_flow
                    if ix_f == ix_t
                        J[F_ix_f_r, X_ix_t_fst] = -1.0
                        J[F_ix_f_i, X_ix_t_snd] = -1.0
                    end
                elseif b == PSY.BusTypes.PV
                    # State variables are Reactive Power Generated and Voltage Angle
                    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - p_gen[i]
                    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
                    # x[2*i] (associated with q_gen) does not appear in q_flow
                    # x[2*i] (associated with q_gen) does not appear in the active power balance
                    if ix_f == ix_t
                        #Jac: Reactive PF against local active power
                        J[F_ix_f_i, X_ix_t_fst] = -1.0
                        #Jac: Active PF against same Angle: θ[ix_f] =  θ[ix_t]
                        J[F_ix_f_r, X_ix_t_snd] =
                            Vm[ix_f] * sum(
                                Vm[k] * (
                                    real(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k]) +
                                    imag(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k])
                                ) for k in neighbors[ix_f] if k != ix_f
                            )
                        #Jac: Reactive PF against same Angle: θ[ix_f] = θ[ix_t]
                        J[F_ix_f_i, X_ix_t_snd] =
                            Vm[ix_f] * sum(
                                Vm[k] * (
                                    real(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k]) -
                                    imag(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k])
                                ) for k in neighbors[ix_f] if k != ix_f
                            )
                    else
                        g_ij = real(Yb[ix_f, ix_t])
                        b_ij = imag(Yb[ix_f, ix_t])
                        #Jac: Active PF against other angles θ[ix_t]
                        J[F_ix_f_r, X_ix_t_snd] =
                            Vm[ix_f] *
                            Vm[ix_t] *
                            (g_ij * sin(θ[ix_f] - θ[ix_t]) + b_ij * -cos(θ[ix_f] - θ[ix_t]))
                        #Jac: Reactive PF against other angles θ[ix_t]
                        J[F_ix_f_i, X_ix_t_snd] =
                            Vm[ix_f] *
                            Vm[ix_t] *
                            (g_ij * -cos(θ[ix_f] - θ[ix_t]) - b_ij * sin(θ[ix_f] - θ[ix_t]))
                    end
                elseif b == PSY.BusTypes.PQ
                    # State variables are Voltage Magnitude and Voltage Angle
                    # Everything appears in everything
                    if ix_f == ix_t
                        #Jac: Active PF against same voltage magnitude Vm[ix_f]
                        J[F_ix_f_r, X_ix_t_fst] =
                            2 * real(Yb[ix_f, ix_t]) * Vm[ix_f] + sum(
                                Vm[k] * (
                                    real(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k]) +
                                    imag(Yb[ix_f, k]) * sin(θ[ix_f] - θ[k])
                                ) for k in neighbors[ix_f] if k != ix_f
                            )
                        #Jac: Active PF against same angle θ[ix_f]
                        J[F_ix_f_r, X_ix_t_snd] =
                            Vm[ix_f] * sum(
                                Vm[k] * (
                                    real(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k]) +
                                    imag(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k])
                                ) for k in neighbors[ix_f] if k != ix_f
                            )

                        #Jac: Reactive PF against same voltage magnitude Vm[ix_f]
                        J[F_ix_f_i, X_ix_t_fst] =
                            -2 * imag(Yb[ix_f, ix_t]) * Vm[ix_f] + sum(
                                Vm[k] * (
                                    real(Yb[ix_f, k]) * sin(θ[ix_f] - θ[k]) -
                                    imag(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k])
                                ) for k in neighbors[ix_f] if k != ix_f
                            )
                        #Jac: Reactive PF against same angle θ[ix_f]
                        J[F_ix_f_i, X_ix_t_snd] =
                            Vm[ix_f] * sum(
                                Vm[k] * (
                                    real(Yb[ix_f, k]) * cos(θ[ix_f] - θ[k]) -
                                    imag(Yb[ix_f, k]) * -sin(θ[ix_f] - θ[k])
                                ) for k in neighbors[ix_f] if k != ix_f
                            )
                    else
                        g_ij = real(Yb[ix_f, ix_t])
                        b_ij = imag(Yb[ix_f, ix_t])
                        #Jac: Active PF w/r to different voltage magnitude Vm[ix_t]
                        J[F_ix_f_r, X_ix_t_fst] =
                            Vm[ix_f] *
                            (g_ij * cos(θ[ix_f] - θ[ix_t]) + b_ij * sin(θ[ix_f] - θ[ix_t]))
                        #Jac: Active PF w/r to different angle θ[ix_t]
                        J[F_ix_f_r, X_ix_t_snd] =
                            Vm[ix_f] *
                            Vm[ix_t] *
                            (g_ij * sin(θ[ix_f] - θ[ix_t]) + b_ij * -cos(θ[ix_f] - θ[ix_t]))

                        #Jac: Reactive PF w/r to different voltage magnitude Vm[ix_t]
                        J[F_ix_f_i, X_ix_t_fst] =
                            Vm[ix_f] *
                            (g_ij * sin(θ[ix_f] - θ[ix_t]) - b_ij * cos(θ[ix_f] - θ[ix_t]))
                        #Jac: Reactive PF w/r to different angle θ[ix_t]
                        J[F_ix_f_i, X_ix_t_snd] =
                            Vm[ix_f] *
                            Vm[ix_t] *
                            (g_ij * -cos(θ[ix_f] - θ[ix_t]) - b_ij * sin(θ[ix_f] - θ[ix_t]))
                    end
                else
                    error("Undefined Conditional")
                end
            end
        end
    end

    res = similar(x0)
    pf!(res, x0)
    if sum(res) > 10 * (N_BUS * 2)
        _, ix = findmax(res)
        bus_no = PSY.get_number(buses[(ix + 1) ÷ 2])
        @warn "Initial guess provided results in a large initial residual. Largest residual at bus $bus_no"
    end

    if finite_diff
        res = NLsolve.nlsolve(pf!, x0; kwargs...)
    else
        F0 = similar(x0)
        df = NLsolve.OnceDifferentiable(pf!, jsp!, x0, F0, J0)
        res = NLsolve.nlsolve(df, x0; kwargs...)
    end

    return res
end
