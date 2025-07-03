const PTDFPowerFlowData = PowerFlowData{
    PNM.PTDF{
        Tuple{Vector{Int64}, Vector{String}},
        Tuple{Dict{Int64, Int64}, Dict{String, Int64}},
        Matrix{Float64},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

const vPTDFPowerFlowData = PowerFlowData{
    PNM.VirtualPTDF{
        Tuple{Vector{String}, Vector{Int64}},
        Tuple{Dict{String, Int64}, Dict{Int64, Int64}},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

const ABAPowerFlowData = PowerFlowData{
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
    PNM.BA_Matrix{
        Tuple{Vector{Int64}, Vector{String}},
        Tuple{Dict{Int64, Int64}, Dict{String, Int64}}},
}

"""
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::PTDFPowerFlowData`:
        PTDFPowerFlowData structure containing all the information related to the system's power flow.
"""
function solve_power_flow_data!(
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
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::vPTDFPowerFlowData`:
        vPTDFPowerFlowData structure containing all the information related to the system's power flow.
"""
function solve_power_flow_data!(
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
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::ABAPowerFlowData`:
        ABAPowerFlowData structure containing all the information related to the system's power flow.
"""
# DC flow: ABA and BA case
function solve_power_flow_data!(
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
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
- `::PTDFDCPowerFlow`:
        use PTDFDCPowerFlow() to evaluate the power flows according to the
        method based on the PTDF matrix
"""
function solve_power_flow(
    sys::PSY.System,
    ::PTDFDCPowerFlow;
    correct_bustypes = false,
)
    data = PowerFlowData(PTDFDCPowerFlow(), sys; correct_bustypes = correct_bustypes)
    solve_power_flow_data!(data)
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
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
- `::DCPowerFlow`:
        use DCPowerFlow() to evaluate the power flows according to the method
        based on the ABA and BA matrices
"""
function solve_power_flow(
    sys::PSY.System,
    ::DCPowerFlow;
    correct_bustypes = false,
)
    data = PowerFlowData(DCPowerFlow(), sys; correct_bustypes = correct_bustypes)
    solve_power_flow_data!(data)
    return write_results(data, sys)
end

"""
Evaluates the power flows on each system's branch by means of the Virtual PTDF
matrix.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the number of timestep considered in "data".
The DataFrame containts the flows and angles related to the information stored
in the PSY.System considered as input.

# Arguments:
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
- `::vPTDFDCPowerFlow`:
        use vPTDFDCPowerFlow() to evaluate the power flows according to the
        method based on the Virtual PTDF matrix
"""
function solve_power_flow(
    sys::PSY.System,
    ::vPTDFDCPowerFlow;
    correct_bustypes = false,
)
    data = PowerFlowData(vPTDFDCPowerFlow(), sys; correct_bustypes = correct_bustypes)
    solve_power_flow_data!(data)
    return write_results(data, sys)
end

# MULTI PERIOD ###############################################################

"""
Evaluates the power flows on each system's branch by means of the PTDF matrix.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the number of timestep considered in "data".
Each DataFrame containts the flows and angles.

# Arguments:
- `data::PTDFPowerFlowData`:
        PowerFlowData structure containing the system's data per each timestep
        considered, as well as the PTDF matrix.
- `sys::PSY.System`:
        container gathering the system data.
"""
function solve_power_flow_data(
    data::PTDFPowerFlowData,
    sys::PSY.System;
)
    solve_power_flow_data!(data)
    return write_results(data, sys)
end

"""
Evaluates the power flows on each system's branch by means of the ABA and BA
matrices.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the number of timestep considered in "data".
Each DataFrame containts the flows and angles.

# Arguments:
- `data::ABAPowerFlowData`:
        PowerFlowData structure containing the system's data per each timestep
        considered, as well as the ABA and BA matrices.
- `sys::PSY.System`:
        container gathering the system data.
"""
function solve_power_flow_data(
    data::ABAPowerFlowData,
    sys::PSY.System;
)
    solve_power_flow_data!(data)
    return write_results(data, sys)
end

"""
Evaluates the power flows on each system's branch by means of Virtual PTDF
matrices.
Updates the PowerFlowData structure "data" and returns a dictionary containing
a number of DataFrames equal to the number of timestep considered in "data".
Each DataFrame containts the flows and angles.

# Arguments:
- `data::PTDFPowerFlowData`:
        PowerFlowData structure containing the system data per each timestep
        considered, as well as the Virtual PTDF matrix.
- `sys::PSY.System`:
        container gathering the system data.
"""
function solve_power_flow_data(
    data::vPTDFPowerFlowData,
    sys::PSY.System;
)
    solve_power_flow_data!(data)
    return write_results(data, sys)
end
