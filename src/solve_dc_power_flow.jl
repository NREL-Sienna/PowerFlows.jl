"""
Adjust the power injections vector to account for the power flows through LCCs.
    
Relies on the fact that we calculate those flows during initialization and save them
to the `active_power_flow_from_to` and `active_power_flow_to_from` fields of the
`LCCParameters` struct.
"""
function adjust_power_injections_for_lccs!(power_injections::Matrix{Float64},
    lcc_params::LCCParameters,
)
    for (i, bus_inds) in enumerate(lcc_params.bus_indices)
        from_bus_ix, to_bus_ix = bus_inds
        rectifier_power = lcc_params.arc_active_power_flow_from_to[i]
        # inverter_power here takes into account losses.
        inverter_power = lcc_params.arc_active_power_flow_to_from[i]
        power_injections[from_bus_ix, :] .-= rectifier_power
        power_injections[to_bus_ix, :] .+= inverter_power
    end
    return
end

"""
    solve_power_flow!(data::PTDFPowerFlowData)
Evaluates the PTDF power flow and writes the result to the fields of the 
[`PTDFPowerFlowData`](@ref) structure.

This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_active_power_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_active_power_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
function solve_power_flow!(
    data::PTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    # get net power injections
    power_injections = data.bus_active_power_injections .- data.bus_active_power_withdrawals
    power_injections .+= data.bus_hvdc_net_power
    # evaluate flows
    data.arc_active_power_flow_from_to .=
        data.power_network_matrix.data' * power_injections
    data.arc_active_power_flow_to_from .= -data.arc_active_power_flow_from_to
    # HVDC flows stored separately and already calculated: see initialize_power_flow_data!
    valid_ix = get_valid_ix(data)
    p_inj = power_injections[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    _compute_branch_angle_differences_from_data!(data)
    data.converged .= true
    if get_calculate_loss_factors(data)
        data.loss_factors .= dc_loss_factors(data, power_injections)
    end
    return
end

"""
    solve_power_flow!(data::vPTDFPowerFlowData)

Evaluates the virtual PTDF power flow and writes the results to the fields 
of the [`vPTDFPowerFlowData`](@ref) structure.


This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_active_power_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_active_power_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
function solve_power_flow!(
    data::vPTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    power_injections = data.bus_active_power_injections .- data.bus_active_power_withdrawals
    power_injections .+= data.bus_hvdc_net_power
    data.arc_active_power_flow_from_to .=
        my_mul_mt(data.power_network_matrix, power_injections)
    data.arc_active_power_flow_to_from .= -data.arc_active_power_flow_from_to
    # HVDC flows stored separately and already calculated: see initialize_power_flow_data!
    valid_ix = get_valid_ix(data)
    p_inj = power_injections[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    _compute_branch_angle_differences_from_data!(data)
    data.converged .= true
    if get_calculate_loss_factors(data)
        data.loss_factors .= dc_loss_factors(data, power_injections)
    end
    return
end

# TODO: solve just for some lines with vPTDF

"""
    solve_power_flow!(data::ABAPowerFlowData)

Evaluates the DC power flow and writes the results (branch flows) to the fields 
of the [`ABAPowerFlowData`](@ref) structure.


This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_active_power_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_active_power_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
# DC flow: ABA and BA case
function solve_power_flow!(
    data::ABAPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.power_network_matrix.data)
    full_factor!(solver_cache, data.power_network_matrix.data)
    # get net injections
    power_injections = data.bus_active_power_injections - data.bus_active_power_withdrawals
    power_injections .+= data.bus_hvdc_net_power
    # save angles and power flows
    valid_ix = get_valid_ix(data)
    p_inj = power_injections[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    data.arc_active_power_flow_from_to .= data.aux_network_matrix.data' * data.bus_angles
    data.arc_active_power_flow_to_from .= -data.arc_active_power_flow_from_to
    # HVDC flows stored separately and already calculated: see initialize_power_flow_data!
    _compute_branch_angle_differences_from_data!(data)
    data.converged .= true
    return
end

# SINGLE PERIOD ##############################################################

"""
    solve_power_flow(
        pf::T,
        sys::PSY.System,
        flow_reporting::FlowReporting
    ) where T <: AbstractDCPowerFlow


Evaluates the provided DC power flow method `pf` on the [PowerSystems.System](@extref) `sys`,
returning a dictionary of `DataFrame`s containing the calculated flows and bus angles.
The flow_reporting input determines if flows are reported for arcs (FlowReporting.ARC_FLOWS)
or for branches (FlowReporting.BRANCH_FLOWS)

Configuration options like `time_steps`, `time_step_names`, `network_reductions`, and
`correct_bustypes` should be set on the power flow object (e.g., `DCPowerFlow(; time_steps=2)`).

Provided for convenience: this interface bypasses the need to create a `PowerFlowData`
struct, but that's still what's happening under the hood.

# Example
```julia
using PowerFlows, PowerSystemCaseBuilder
sys = build_system(PSITestSystems, "c_sys5")
d = solve_power_flow(DCPowerFlow(), sys)
display(d["1"]["flow_results"])
display(d["1"]["bus_results"])
```
"""
function solve_power_flow(
    pf::T,
    sys::PSY.System,
    flow_reporting::FlowReporting,
) where {T <: AbstractDCPowerFlow}
    with_units_base(sys, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(pf, sys)
        solve_power_flow!(data)
        return write_results(data, sys, flow_reporting)
    end
end

# MULTI PERIOD ###############################################################

"""
Evaluates the power flows on the system's branches by means of the method associated with
the `PowerFlowData` structure `data`, which can be one of `PTDFPowerFlowData`,
`vPTDFPowerFlowData`, or `ABAPowerFlowData`.
Returns a dictionary of `DataFrame`s, each containing the flows and bus voltages for
the input `PSY.System` at that time_step.
The flow_reporting input determines if flows are reported for arcs (FlowReporting.ARC_FLOWS)
or for branches (FlowReporting.BRANCH_FLOWS)

# Arguments:
- `data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData}`:
        `PowerFlowData` structure containing the system's data per each time_step
        considered, as well as the associated matrix for the power flow.
- `sys::PSY.System`:
        container gathering the system data.
- `flow_reporting::FlowReporting`:
        Format for reporting flows

Note that `data` must have been created from the [System](@extref PowerSystems.System) 
`sys` using one of the [`PowerFlowData`](@ref) constructors.

# Example
```julia
using PowerFlows, PowerSystemCaseBuilder
sys = build_system(PSITestSystems, "c_sys14")
data = PowerFlowData(PTDFDCPowerFlow(; time_steps = 2), sys)
d = solve_power_flow(data, sys)
display(d["2"]["flow_results"])
```
"""
function solve_power_flow(
    data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData},
    sys::PSY.System,
    flow_reporting::FlowReporting;
)
    solve_power_flow!(data)
    return write_results(data, sys, flow_reporting)
end

"""
    _get_arc_resistances(data::Union{PTDFPowerFlowData, vPTDFPowerFlowData}) -> Vector{Float64}

Look up the resistance of each arc from the network reduction data.
"""
function _get_arc_resistances(
    data::Union{PTDFPowerFlowData, vPTDFPowerFlowData},
)
    nrd = get_network_reduction_data(data)
    arc_ax = get_arc_axis(data)
    Rs = zeros(length(arc_ax))
    # TODO simpler way? Should be a uniform interface for this type of thing...
    for (ix_arc, arc) in enumerate(arc_ax)
        if arc in keys(PNM.get_direct_branch_map(nrd))
            line = PNM.get_direct_branch_map(nrd)[arc]
            r = PSY.get_r(line)
        elseif arc in keys(PNM.get_parallel_branch_map(nrd))
            parallel_lines = PNM.parallel_branch_map(nrd)[arc]
            r = 1 / (sum(1 / PSY.get_r.(parallel_lines.branches)))
        elseif arc in keys(PNM.get_series_branch_map(nrd))
            series_lines = PNM.series_branch_map(nrd)[arc]
            r = sum(PSY.get_r.(series_chain) for series_chain in series_lines)
        elseif arc in keys(PNM.get_transformer3W_map(nrd))
            transformer3w = PNM.transformer3W_map(nrd)[arc]
            r = PSY.get_equivalent_r(transformer3w)
        else
            error("Arc $arc not found in any of the branch maps.")
        end
        Rs[ix_arc] = r
    end
    return Rs
end

"""
    dc_loss_factors(
        data::Union{PTDFPowerFlowData, vPTDFPowerFlowData},
        P::Matrix{Float64},
    ) -> Matrix{Float64}

Compute the gradient of total system active power losses with respect to
bus injections using the DC power flow approximation:

    ∂Loss/∂P = 2 · PTDFᵀ · diag(R) · PTDF · P

This is equivalent to the per-element form:

    ∂Loss/∂Pᵢ = Σₖ 2·Rₖ·PTDFₖᵢ·Σⱼ PTDFₖⱼ·Pⱼ

# Arguments
- `data::Union{PTDFPowerFlowData, vPTDFPowerFlowData}`: solved power flow data containing
  the PTDF matrix and network reduction data for looking up branch resistances.
- `P::Matrix{Float64}`: bus injection matrix of size `(num_buses, num_timesteps)`.

# Returns
- `Matrix{Float64}`: loss factor matrix of size `(num_buses, num_timesteps)`, where each
  entry `[i, t]` is the marginal change in total system losses per unit injection at bus `i`
  in time step `t`.
"""
function dc_loss_factors(
    data::PTDFPowerFlowData,
    P::Matrix{Float64},
)
    Rs = _get_arc_resistances(data)
    ptdf_t = data.power_network_matrix.data
    # PERF could be optimized: remove the Diagonal call.
    return 2 * ptdf_t * LinearAlgebra.Diagonal(Rs) * ptdf_t' * P
end

function dc_loss_factors(
    data::vPTDFPowerFlowData,
    P::Matrix{Float64},
)
    Rs = _get_arc_resistances(data)
    ptdf = data.power_network_matrix
    # PTDF * P via row-by-row access (VirtualPTDF has no .data field)
    flows = my_mul_mt(ptdf, P)  # (n_arcs × n_ts)
    weighted_flows = Rs .* flows
    # PTDFᵀ * weighted_flows, accumulated row-by-row
    arc_ax = get_arc_axis(data)
    n_buses = size(P, 1)
    n_ts = size(P, 2)
    result = zeros(n_buses, n_ts)
    for (k, arc) in enumerate(arc_ax)
        row_k = ptdf[arc, :]  # length n_buses
        for t in 1:n_ts
            w = weighted_flows[k, t]
            for j in 1:n_buses
                result[j, t] += row_k[j] * w
            end
        end
    end
    return 2 .* result
end
