"""
    struct ACPowerFlowResidual

A struct to keep track of the residuals in the Newton-Raphson AC power flow calculation.

# Fields
- `data::ACPowerFlowData`: The grid model data.
- `Rf!::Function`: A function that updates the residuals based on the latest values stored in the grid at the given iteration.
- `Rv::Vector{Float64}`: A vector of the values of the residuals.
- `P_net::Vector{Float64}`: A vector of net active power injections.
- `Q_net::Vector{Float64}`: A vector of net reactive power injections.
- `P_net_set::Vector{Float64}`: A vector of the set-points for active power injections (their initial values before power flow calculation).
- `bus_slack_participation_factors::SparseVector{Float64, Int}`: A sparse vector of the slack participation factors aggregated at the bus level.
- `subnetworks::Dict{Int64, Vector{Int64}}`: The dictionary that identifies subnetworks (connected components), with the key defining the REF bus, values defining the corresponding buses in the subnetwork.
"""
struct ACPowerFlowResidual
    data::ACPowerFlowData
    Rf!::Function
    Rv::Vector{Float64}
    P_net::Vector{Float64}
    Q_net::Vector{Float64}
    P_net_set::Vector{Float64}
    bus_slack_participation_factors::SparseVector{Float64, Int}
    subnetworks::Dict{Int64, Vector{Int64}}
end

"""
    ACPowerFlowResidual(data::ACPowerFlowData, time_step::Int64)

Create an instance of `ACPowerFlowResidual` for a given time step.

# Arguments
- `data::ACPowerFlowData`: The power flow data representing the power system model.
- `time_step::Int64`: The time step for which the power flow calculation is executed.

# Returns
- `ACPowerFlowResidual`: An instance containing the residual values, net bus active power injections,
and net bus reactive power injections.
"""
function ACPowerFlowResidual(data::ACPowerFlowData, time_step::Int64)
    n_buses = first(size(data.bus_type))
    n_lccs = size(data.lcc.p_set, 1)
    P_net = Vector{Float64}(undef, n_buses)
    Q_net = Vector{Float64}(undef, n_buses)

    P_net_set = zeros(Float64, n_buses)
    bus_type = view(data.bus_type, :, time_step)

    spf_idx = Int[]
    spf_val = Float64[]
    sum_sl_weights = 0.0  # for scope

    # ref_bus is set to the first REF bus found - will be used for the total slack power
    subnetworks =
        _find_subnetworks_for_reference_buses(data.power_network_matrix.data, bus_type)

    for (ix, bt) in zip(1:n_buses, bus_type)
        P_net[ix] =
            data.bus_activepower_injection[ix, time_step] -
            get_bus_activepower_total_withdrawals(data, ix, time_step)
        Q_net[ix] =
            data.bus_reactivepower_injection[ix, time_step] -
            get_bus_reactivepower_total_withdrawals(data, ix, time_step)
        P_net_set[ix] = P_net[ix]

        bt ∈ (PSY.ACBusTypes.REF, PSY.ACBusTypes.PV) || continue
        (spf_v = data.bus_slack_participation_factors[ix, time_step]) == 0.0 && continue
        push!(spf_idx, ix)
        push!(spf_val, spf_v)
        sum_sl_weights += spf_v
    end

    if sum_sl_weights == 0.0
        throw(ArgumentError("sum of slack_participation_factors cannot be zero"))
    end

    if any(spf_val .< 0.0)
        throw(ArgumentError("slack_participation_factors cannot be negative"))
    end

    # Actually should be fine e.g. when PV is changed to PQ
    # if sum_sl_weights != sum(abs.(data.slack_participation_factors[:, time_step]))
    #     warn("Only slack weights for REF and PV buses can be considered.")
    # end

    # bus slack participation factors relevant for the current time step:
    bus_slack_participation_factors = sparsevec(spf_idx, spf_val, n_buses)

    # normalize slack participation factors to sum to 1 per every subnetwork
    for subnetwork_buses in values(subnetworks)
        bspf_subnetwork = view(bus_slack_participation_factors, subnetwork_buses)
        sum_bspf_subnetwork = sum(bspf_subnetwork)
        sum_bspf_subnetwork == 0.0 && throw(
            ArgumentError(
                "sum of slack_participation_factors per subnetwork cannot be zero",
            ),
        )
        bspf_subnetwork ./= sum_bspf_subnetwork
    end

    return ACPowerFlowResidual(
        data,
        _update_residual_values!,
        Vector{Float64}(undef, 2 * n_buses + 4 * n_lccs),
        P_net,
        Q_net,
        P_net_set,
        bus_slack_participation_factors,
        subnetworks,
    )
end

"""
    (Residual::ACPowerFlowResidual)(Rv::Vector{Float64}, x::Vector{Float64}, time_step::Int64)

Evaluate the AC power flow residuals and store the result in `Rv` using the provided
state vector `x` and the current time step `time_step`.
The residuals are updated inplace in the struct and additionally copied to the provided array.
This function implements the functor approach for the `ACPowerFlowResidual` struct.
This makes the struct callable.
Calling the `ACPowerFlowResidual` will also update the values of P, Q, V, Θ in the `data` struct.

# Arguments
- `Rv::Vector{Float64}`: The vector to store the calculated residuals.
- `x::Vector{Float64}`: The state vector.
- `time_step::Int64`: The current time step.
"""
function (Residual::ACPowerFlowResidual)(
    Rv::Vector{Float64},
    x::Vector{Float64},
    time_step::Int64,
)
    Residual.Rf!(
        Residual.Rv,
        x,
        Residual.P_net,
        Residual.Q_net,
        Residual.P_net_set,
        Residual.bus_slack_participation_factors,
        Residual.subnetworks,
        Residual.data,
        time_step,
    )
    copyto!(Rv, Residual.Rv)
    return
end

"""
    (Residual::ACPowerFlowResidual)(x::Vector{Float64}, time_step::Int64)

Update the AC power flow residuals inplace and store the result in the attribute `Rv` of the struct.
The inputs are the values of state vector `x` and the current time step `time_step`.
This function implements the functor approach for the `ACPowerFlowResidual` struct.
This makes the struct callable.
Calling the `ACPowerFlowResidual` will also update the values of P, Q, V, Θ in the `data` struct.

# Arguments
- `x::Vector{Float64}`: The state vector values.
- `time_step::Int64`: The current time step.
"""
function (Residual::ACPowerFlowResidual)(x::Vector{Float64}, time_step::Int64)
    Residual.Rf!(
        Residual.Rv,
        x,
        Residual.P_net,
        Residual.Q_net,
        Residual.P_net_set,
        Residual.bus_slack_participation_factors,
        Residual.subnetworks,
        Residual.data,
        time_step,
    )
    return
end

function _setpq(
    ix::Int,
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
)
    # Set the active and reactive power injections at the bus
    data.bus_activepower_injection[ix, time_step] =
        P_net[ix] + get_bus_activepower_total_withdrawals(data, ix, time_step)
    data.bus_reactivepower_injection[ix, time_step] =
        Q_net[ix] + get_bus_reactivepower_total_withdrawals(data, ix, time_step)
end

# dispatching on Val for performance reasons.
function _set_state_variables_at_bus!(
    ix::Int,
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    P_net_set::Vector{Float64},
    P_slack::Float64,
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.REF})
    # When bustype == REFERENCE PSY.ACACBus, state variables are Active and Reactive Power Generated
    P_net[ix] = P_net_set[ix] + P_slack
    Q_net[ix] = StateVector[2 * ix]
    _setpq(
        ix,
        P_net,
        Q_net,
        data,
        time_step,
    )
end

function _set_state_variables_at_bus!(
    ix::Int,
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    P_net_set::Vector{Float64},
    P_slack::Float64,
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.PV})
    # When bustype == PV PSY.ACACBus, state variables are Reactive Power Generated and Voltage Angle
    # We still update both P and Q values in case the PV bus participates in distributed slack
    P_net[ix] = P_net_set[ix] + P_slack
    Q_net[ix] = StateVector[2 * ix - 1]
    _setpq(
        ix,
        P_net,
        Q_net,
        data,
        time_step,
    )
    data.bus_angles[ix, time_step] = StateVector[2 * ix]
end

function _set_state_variables_at_bus!(
    ix::Int,
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    ::Vector{Float64},
    ::Float64,
    StateVector::Vector{Float64},
    data::ACPowerFlowData,
    time_step::Int64,
    ::Val{PSY.ACBusTypes.PQ})
    # When bustype == PQ PSY.ACBus, state variables are Voltage Magnitude and Voltage Angle
    # delta_vm = (vm_1 = data.bus_magnitude[ix, time_step]) - StateVector[2 * ix - 1]
    vm_1 = data.bus_magnitude[ix, time_step]
    vm_2 = StateVector[2 * ix - 1]
    data.bus_magnitude[ix, time_step] = vm_2
    data.bus_angles[ix, time_step] = StateVector[2 * ix]
    # update P_net and Q_net for ZIP loads
    P_net[ix] +=
        data.bus_activepower_constant_current_withdrawals[ix, time_step] * (vm_1 - vm_2) +
        data.bus_activepower_constant_impedance_withdrawals[ix, time_step] *
        (vm_1^2 - vm_2^2)
    Q_net[ix] +=
        data.bus_reactivepower_constant_current_withdrawals[ix, time_step] * (vm_1 - vm_2) +
        data.bus_reactivepower_constant_impedance_withdrawals[ix, time_step] *
        (vm_1^2 - vm_2^2)
    # set the active and reactive power injections at the bus
    _setpq(
        ix,
        P_net,
        Q_net,
        data,
        time_step,
    )
end

"""
    _update_residual_values!(
        F::Vector{Float64},
        x::Vector{Float64},
        P_net::Vector{Float64},
        Q_net::Vector{Float64},
        data::ACPowerFlowData,
        time_step::Int64,
    )

Update the residual values for the Newton-Raphson AC power flow calculation. This function is used internally in the
`ACPowerFlowResidual` struct. This function also updates the values of P, Q, V, Θ in the `data` struct.

# Arguments
- `F::Vector{Float64}`: Vector of the values of the residuals.
- `x::Vector{Float64}`: State vector values.
- `P_net::Vector{Float64}`: Vector of net active power injections at each bus.
- `Q_net::Vector{Float64}`: Vector of net reactive power injections at each bus.
- `P_net_set::Vector{Float64}`: Vector of the set-points for active power injections (their initial values before power flow calculation).
- `bus_slack_participation_factors::SparseVector{Float64, Int}`: Sparse vector of the slack participation factors aggregated at the bus level.
- `ref_bus::Int`: The index of the reference bus to be used for the total slack power.
- `data::ACPowerFlowData`: Data structure representing the grid model for the AC power flow calculation.
- `time_step::Int64`: The current time step for which the residual values are being updated.
"""
function _update_residual_values!(
    F::Vector{Float64},
    x::Vector{Float64},
    P_net::Vector{Float64},
    Q_net::Vector{Float64},
    P_net_set::Vector{Float64},
    bus_slack_participation_factors::SparseVector{Float64, Int},
    subnetworks::Dict{Int64, Vector{Int64}},
    data::ACPowerFlowData,
    time_step::Int64,
)
    # update P_net, Q_net, data.bus_angles, data.bus_magnitude based on X
    Yb = data.power_network_matrix.data
    Yb_facts = data.power_network_matrix.data_facts
    num_lcc = size(data.lcc.p_set, 1)
    bus_types = view(data.bus_type, :, time_step)

    for (ref_bus, subnetwork_buses) in subnetworks
        P_slack =
            (x[2 * ref_bus - 1] - P_net_set[ref_bus]) .*
            bus_slack_participation_factors[subnetwork_buses]

        for (ix, bt, p_bus_slack) in
            zip(subnetwork_buses, bus_types[subnetwork_buses], P_slack)
            # creating Val(bt) at runtime is slow, requires allocating: split into cases
            # explicitly, so instead it's Val(compile-time constant).
            if bt == PSY.ACBusTypes.PQ
                _set_state_variables_at_bus!(
                    ix,
                    P_net,
                    Q_net,
                    P_net_set,
                    p_bus_slack,
                    x,
                    data,
                    time_step,
                    Val(PSY.ACBusTypes.PQ),
                )
            elseif bt == PSY.ACBusTypes.PV
                _set_state_variables_at_bus!(
                    ix,
                    P_net,
                    Q_net,
                    P_net_set,
                    p_bus_slack,
                    x,
                    data,
                    time_step,
                    Val(PSY.ACBusTypes.PV),
                )
            elseif bt == PSY.ACBusTypes.REF
                _set_state_variables_at_bus!(
                    ix,
                    P_net,
                    Q_net,
                    P_net_set,
                    p_bus_slack,
                    x,
                    data,
                    time_step,
                    Val(PSY.ACBusTypes.REF),
                )
            end
        end
    end

    if num_lcc > 0
        data.lcc.rectifier_tap[:, time_step] = x[(end - 4 * num_lcc + 1):4:end]
        data.lcc.inverter_tap[:, time_step] = x[(end - 4 * num_lcc + 2):4:end]
        data.lcc.rectifier_delay_angle[:, time_step] = x[(end - 4 * num_lcc + 3):4:end]
        data.lcc.inverter_extinction_angle[:, time_step] = x[(end - 4 * num_lcc + 4):4:end]
        _update_ybus_lcc!(Yb_facts, data, time_step)
    end

    # compute active, reactive power balances using the just updated values.
    Vm = view(data.bus_magnitude, :, time_step)
    θ = view(data.bus_angles, :, time_step)
    # F is active and reactive power balance equations at all buses
    Yb_vals = SparseArrays.nonzeros(Yb)
    Yb_rowvals = SparseArrays.rowvals(Yb)
    F .= 0.0
    for bus_to in axes(Yb, 1)
        for j in Yb.colptr[bus_to]:(Yb.colptr[bus_to + 1] - 1)
            yb = Yb_vals[j]
            bus_from = Yb_rowvals[j]
            gb = real(yb)
            bb = imag(yb)
            Δθ = θ[bus_from] - θ[bus_to]
            if bus_from == bus_to
                F[2 * bus_from - 1] += Vm[bus_from] * Vm[bus_to] * gb
                F[2 * bus_from] += -Vm[bus_from] * Vm[bus_to] * bb
            else
                F[2 * bus_from - 1] +=
                    Vm[bus_from] * Vm[bus_to] * (gb * cos(Δθ) + bb * sin(Δθ))
                F[2 * bus_from] += Vm[bus_from] * Vm[bus_to] * (gb * sin(Δθ) - bb * cos(Δθ))
            end
        end
    end

    F[1:2:(end - 4 * num_lcc)] .-= P_net
    F[2:2:(end - 4 * num_lcc)] .-= Q_net

    if num_lcc > 0
        P_lcc_from =
            Vm[data.lcc.rectifier_bus, time_step] .* data.lcc.rectifier_tap[:, time_step] .*
            sqrt(6) / π .* data.lcc.rectifier_i_dc[:, time_step] .*
            cos.(data.lcc.rectifier_delay_angle[:, time_step]) .-
            sqrt(3 / 2) * sqrt(6) / π .* data.lcc.rectifier_transformer_reactance .*
            data.lcc.rectifier_i_dc[:, time_step] .^ 2
        P_lcc_to =
            Vm[data.lcc.inverter_bus, time_step] .* data.lcc.inverter_tap[:, time_step] .*
            sqrt(6) / π .* data.lcc.inverter_i_dc[:, time_step] .*
            cos.(data.lcc.inverter_extinction_angle[:, time_step]) .-
            sqrt(3 / 2) * sqrt(6) / π .* data.lcc.inverter_transformer_reactance .*
            data.lcc.inverter_i_dc[:, time_step] .^ 2
        F[(end - 4 * num_lcc + 1):4:end] .= P_lcc_from .- data.lcc.p_set[:, time_step]
        F[(end - 4 * num_lcc + 2):4:end] .=
            P_lcc_from .+ P_lcc_to .-
            data.lcc.dc_line_resistance .* data.lcc.rectifier_i_dc[:, time_step] .^ 2
        F[(end - 4 * num_lcc + 3):4:end] .=
            data.lcc.rectifier_delay_angle[:, time_step] .- data.lcc.rectifier_min_alpha
        F[(end - 4 * num_lcc + 4):4:end] .=
            data.lcc.inverter_extinction_angle[:, time_step] .- data.lcc.inverter_min_gamma
    end
    return
end

function _find_subnetworks_for_reference_buses(
    Ybus::SparseMatrixCSC,
    bus_type::AbstractArray{PSY.ACBusTypes},
)
    subnetworks = PNM.find_subnetworks(Ybus, collect(eachindex(bus_type)))
    ref_buses = findall(x -> x == PSY.ACBusTypes.REF, bus_type)
    bus_groups = Dict{Int, Vector{Int}}()
    for (bus_key, subnetwork_buses) in subnetworks
        ref_bus = intersect(ref_buses, subnetwork_buses)
        if length(ref_bus) >= 1
            bus_groups[first(ref_bus)] = collect(subnetwork_buses)
        else
            throw(
                ArgumentError(
                    "No REF bus found in the subnetwork with $(length(subnetwork_buses)) buses defined by bus key $bus_key",
                ),
            )
        end
    end
    return bus_groups
end
