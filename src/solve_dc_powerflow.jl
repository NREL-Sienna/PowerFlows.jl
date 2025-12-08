"""
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::PTDFPowerFlowData`:
        PTDFPowerFlowData structure containing all the information related to the system's power flow.
"""
function solve_powerflow!(
    data::PTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    power_injection .+= data.bus_hvdc_net_power
    # evaluate flows
    data.arc_activepower_flow_from_to .=
        data.power_network_matrix.data' * power_injection
    data.arc_activepower_flow_to_from .= -data.arc_activepower_flow_from_to
    # HVDC flows stored separately and already calculated: see initialize_powerflow_data!
    valid_ix = get_valid_ix(data)
    p_inj = power_injection[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    data.converged .= true
    return
end

"""
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::vPTDFPowerFlowData`:
        vPTDFPowerFlowData structure containing all the information related to the system's power flow.
"""
function solve_powerflow!(
    data::vPTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    power_injection .+= data.bus_hvdc_net_power
    data.arc_activepower_flow_from_to .=
        my_mul_mt(data.power_network_matrix, power_injection)
    data.arc_activepower_flow_to_from .= -data.arc_activepower_flow_from_to
    # HVDC flows stored separately and already calculated: see initialize_powerflow_data!
    valid_ix = get_valid_ix(data)
    p_inj = power_injection[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    data.converged .= true
    return
end

# TODO: solve just for some lines with vPTDF

"""
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::ABAPowerFlowData`:
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
    power_injection .+= data.bus_hvdc_net_power
    # save angles and power flows
    valid_ix = get_valid_ix(data)
    p_inj = power_injection[valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[valid_ix, :] .= p_inj
    data.arc_activepower_flow_from_to .= data.aux_network_matrix.data' * data.bus_angles
    data.arc_activepower_flow_to_from .= -data.arc_activepower_flow_from_to
    # HVDC flows stored separately and already calculated: see initialize_powerflow_data!
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
