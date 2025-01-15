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
function solve_powerflow!(
    data::PTDFPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    # evaluate flows
    data.branch_activepower_flow_from_to .=
        data.power_network_matrix.data' * power_injection
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
    # evaluate bus angles
    p_inj = power_injection[data.valid_ix, :]
    data.bus_angles[data.valid_ix, :] .= data.aux_network_matrix.K \ p_inj
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
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    for i in axes(power_injection, 2)
        # evaluate flows (next line evaluates both both PTDF rows and line flows)
        data.branch_activepower_flow_from_to[:, i] .=
            my_mul_mt(data.power_network_matrix, power_injection[:, i])
        data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
        # evaluate bus angles
        p_inj = power_injection[data.valid_ix, i]
        data.bus_angles[data.valid_ix, i] .= data.aux_network_matrix.K \ p_inj
    end
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
    # get net injections
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    # save angles and power flows
    data.bus_angles[data.valid_ix, :] .=
        data.power_network_matrix.K \ @view power_injection[data.valid_ix, :]
    data.branch_activepower_flow_from_to .= data.aux_network_matrix.data' * data.bus_angles
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
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
- `::PTDFDCPowerFlow`:
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
- `::vPTDFDCPowerFlow`:
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
