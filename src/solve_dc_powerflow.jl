"""
Adjust the power injection vector to account for the power flows through LCCs.
    
Relies on the fact that we calculate those flows during initialization and save them
to the `active_powerflow_from_to` and `active_powerflow_to_from` fields of the
`LCCParameters` struct.
"""
function adjust_power_injection_for_lccs!(power_injection::Matrix{Float64},
    lcc_params::LCCParameters,
)
    for (i, bus_inds) in enumerate(lcc_params.bus_indices)
        from_bus_ix, to_bus_ix = bus_inds
        rectifier_power = lcc_params.arc_activepower_flow_from_to[i]
        # inverter_power here takes into account losses.
        inverter_power = lcc_params.arc_activepower_flow_to_from[i]
        power_injection[from_bus_ix, :] .-= rectifier_power
        power_injection[to_bus_ix, :] .+= inverter_power
    end
    return
end

"""
    solve_powerflow!(data::PTDFPowerFlowData)

Evaluates the PTDF power flow and writes the result to the fields of the 
[`PTDFPowerFlowData`](@ref) structure.

This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_activepower_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_activepower_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
function solve_powerflow!(
    data::PTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    get_lcc_count(data) > 0 && adjust_power_injection_for_lccs!(power_injection, data.lcc)
    # evaluate flows
    data.arc_activepower_flow_from_to .=
        data.power_network_matrix.data' * power_injection
    data.arc_activepower_flow_to_from .= -data.arc_activepower_flow_from_to
    if get_lcc_count(data) > 0
        data.lcc.arc_activepower_flow_to_from .= -data.lcc.arc_activepower_flow_from_to
    end
    # evaluate bus angles
    valid_ix = get_valid_ix(data)
    p_inj = power_injection[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    data.converged .= true
    return
end

"""
    solve_powerflow!(data::vPTDFPowerFlowData)

Evaluates the virtual PTDF power flow and writes the results to the fields 
of the [`vPTDFPowerFlowData`](@ref) structure.


This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_activepower_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_activepower_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
function solve_powerflow!(
    data::vPTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    get_lcc_count(data) > 0 && adjust_power_injection_for_lccs!(power_injection, data.lcc)
    data.arc_activepower_flow_from_to .=
        my_mul_mt(data.power_network_matrix, power_injection)
    data.arc_activepower_flow_to_from .= -data.arc_activepower_flow_from_to
    if get_lcc_count(data) > 0
        data.lcc.arc_activepower_flow_to_from .= -data.lcc.arc_activepower_flow_from_to
    end
    valid_ix = get_valid_ix(data)
    p_inj = power_injection[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    data.converged .= true
    return
end

# TODO: solve just for some lines with vPTDF

"""
    solve_powerflow!(data::ABAPowerFlowData)

Evaluates the DC power flow and writes the results (branch flows) to the fields 
of the [`ABAPowerFlowData`](@ref) structure.


This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_activepower_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_activepower_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
# DC flow: ABA and BA case
function solve_powerflow!(
    data::ABAPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.power_network_matrix.data)
    full_factor!(solver_cache, data.power_network_matrix.data)
    # get net injections
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    get_lcc_count(data) > 0 && adjust_power_injection_for_lccs!(power_injection, data.lcc)
    # save angles and power flows
    valid_ix = get_valid_ix(data)
    p_inj = power_injection[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    data.arc_activepower_flow_from_to .= data.aux_network_matrix.data' * data.bus_angles
    data.arc_activepower_flow_to_from .= -data.arc_activepower_flow_from_to
    if get_lcc_count(data) > 0
        data.lcc.arc_activepower_flow_to_from .= -data.lcc.arc_activepower_flow_from_to
    end
    data.converged .= true
    return
end

# SINGLE PERIOD ##############################################################

"""
    solve_powerflow(
        ::T,
        sys::PSY.System;
    ) where T <: Union{PTDFDCPowerFlow, vPTDFDCPowerFlow, DCPowerFlow}


Evaluates the provided DC power flow method `T` on the `system`, returning a dictionary of 
`DataFrame`s containing the calculated branch flows and bus angles.

Provided for convenience: this interface bypasses the need to create a `PowerFlowData` 
struct, but that's still what's happening under the hood.

# Example
```julia
using PowerFlows, PowerSystemCaseBuilder
sys = PowerSystemCaseBuilder.build_system(PSB.PSITestSystems, "c_sys5")
d = solve_powerflow(DCPowerFlow(), sys)
display(d["1"]["flow_results"])
display(d["1"]["bus_results"])
```
"""
function solve_powerflow(
    ::T,
    sys::PSY.System;
    kargs...,
) where {T <: AbstractDCPowerFlow}
    with_units_base(sys, PSY.UnitSystem.SYSTEM_BASE) do
        data = PowerFlowData(T(), sys; kargs...)
        solve_powerflow!(data)
        return write_results(data, sys)
    end
end

# MULTI PERIOD ###############################################################

"""
Evaluates the power flows on the system's branches by means of the method associated with
the `PowerFlowData` structure `data`, which can be one of `PTDFPowerFlowData`,
`vPTDFPowerFlowData`, or `ABAPowerFlowData`.
Returns a dictionary of `DataFrame`s, each containing the branch flows and bus voltages for
the input `PSY.System` at that timestep.

# Arguments:
- `data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData}`:
        `PowerFlowData` structure containing the system's data per each timestep
        considered, as well as the associated matrix for the power flow.
- `sys::PSY.System`:
        container gathering the system data.

```julia
using PowerFlows, PowerSystemCaseBuilder
sys = PowerSystemCaseBuilder.build_system(PSITestSystems, "c_sys14")
data = PowerFlowData(PTDFDCPowerFlow(), sys, time_steps = 2)
d = solve_powerflow(data, sys)
display(d["2"]["flow_results"])
```
"""
function solve_powerflow(
    data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData},
    sys::PSY.System;
)
    solve_powerflow!(data)
    return write_results(data, sys)
end
