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

Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::PTDFPowerFlowData`:
        [PTDFPowerFlowData](@ref) structure containing all the information related to the system's power flow.
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

Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- [`data::vPTDFPowerFlowData`](@ref vPTDFPowerFlowData):
        a structure containing all the information related to the system's power flow.
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

Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- [`data::ABAPowerFlowData`](@ref ABAPowerFlowData):
        ABAPowerFlowData structure containing all the information related to the system's power flow.
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
Evaluates the power flows on the system's branches by means of the PTDF, virtual PTDF,
or DC power flow method: the type first parameter (a `PTDFDCPowerFlow`, `vPTDFDCPowerFlow`, 
or `DCPowerFlow`) selects the method to be used. Returns a dictionary containing a 
`DataFrame` for the single timestep considered, storing the branch flows and bus 
voltages for the input `PSY.System`.

# Arguments:
- `::Union{PTDFDCPowerFlow, vPTDFDCPowerFlow, DCPowerFlow}`:
        the method of power flow evaluation to be used.
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
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
"""
function solve_powerflow(
    data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData},
    sys::PSY.System;
)
    solve_powerflow!(data)
    return write_results(data, sys)
end
