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
    # evaluate flows
    data.branch_activepower_flow_from_to .=
        data.power_network_matrix.data' * power_injection
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
    # evaluate bus angles
    p_inj = power_injection[data.valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, :] .= p_inj
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
    data.branch_activepower_flow_from_to .=
        my_mul_mt(data.power_network_matrix, power_injection)
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
    p_inj = power_injection[data.valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, :] .= p_inj
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
    # save angles and power flows
    p_inj = power_injection[data.valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, :] .= p_inj
    data.branch_activepower_flow_from_to .= data.aux_network_matrix.data' * data.bus_angles
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
    data.converged .= true
    return
end

# SINGLE PERIOD ##############################################################

"""
Evaluates the power flows on each system's branch by means of the PTDF matrix.
Updates the PowerFlowData structure and returns a dictionary containing a
DataFrame for the single timestep considered.
The DataFrame containts the flows and angles related to the information stored
in the PSY.System considered as input.

# Arguments:
- [`::PTDFDCPowerFlow`](@ref PTDFDCPowerFlow):
        use PTDFDCPowerFlow() to evaluate the power flows according to the
        method based on the PTDF matrix
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
"""
function solve_powerflow(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(PTDFDCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys)
end

"""
Evaluates the power flows on each system's branch by means of the ABA and BA
matrices.
Updates the PowerFlowData structure and returns a dictionary containing a
DataFrame for the single timestep considered.
The DataFrame containts the flows and angles related to the information stored
in the PSY.System considered as input.

# Arguments:
- `::DCPowerFlow`:
        use DCPowerFlow() to evaluate the power flows according to the method
        based on the ABA and BA matrices
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
"""
function solve_powerflow(
    ::DCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(DCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys)
end

"""
Evaluates the power flows on each system's branch by means of the Virtual PTDF
matrix.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the numeber of timestep considered in "data".
The DataFrame containts the flows and angles related to the information stored
in the PSY.System considered as input.

# Arguments:
- [`::vPTDFDCPowerFlow`](@ref vPTDFDCPowerFlow):
        use vPTDFDCPowerFlow() to evaluate the power flows according to the
        method based on the Virtual PTDF matrix
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
"""
function solve_powerflow(
    ::vPTDFDCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(vPTDFDCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys)
end

# MULTI PERIOD ###############################################################

"""
Evaluates the power flows on each system's branch by means of the PTDF matrix.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the numeber of timestep considered in "data".
Each DataFrame containts the flows and angles.

# Arguments:
- `data::PTDFPowerFlowData`:
        PowerFlowData structure containing the system's data per each timestep
        considered, as well as the PTDF matrix.
- `sys::PSY.System`:
        container gathering the system data.
"""
function solve_powerflow(
    data::PTDFPowerFlowData,
    sys::PSY.System;
)
    solve_powerflow!(data)
    return write_results(data, sys)
end

"""
Evaluates the power flows on each system's branch by means of the ABA and BA
matrices.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the numeber of timestep considered in "data".
Each DataFrame containts the flows and angles.

# Arguments:
- `data::ABAPowerFlowData`:
        PowerFlowData structure containing the system's data per each timestep
        considered, as well as the ABA and BA matrices.
- `sys::PSY.System`:
        container gathering the system data.
"""
function solve_powerflow(
    data::ABAPowerFlowData,
    sys::PSY.System;
)
    solve_powerflow!(data)
    return write_results(data, sys)
end

"""
Evaluates the power flows on each system's branch by means of Virtual PTDF
matrices.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the numeber of timestep considered in "data".
Each DataFrame containts the flows and angles.

# Arguments:
- `data::PTDFPowerFlowData`:
        PowerFlowData structure containing the system data per each timestep
        considered, as well as the Virtual PTDF matrix.
- `sys::PSY.System`:
        container gathering the system data.
"""
function solve_powerflow(
    data::vPTDFPowerFlowData,
    sys::PSY.System;
)
    solve_powerflow!(data)
    return write_results(data, sys)
end
