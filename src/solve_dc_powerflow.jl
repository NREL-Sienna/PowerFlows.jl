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

# ? change this to have a more detailed definition ?
const vPTDFPowerFlowData = PowerFlowData{}

const ABAPowerFlowData = PowerFlowData{
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
    PNM.BA_Matrix{
        Tuple{Vector{String}, Vector{Int64}},
        Tuple{Dict{String, Int64}, Dict{Int64, Int64}}},
}

# SINGLE PERIOD: method based on ABA and BA matrices
function _solve_powerflows_single!(
    data::ABAPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    # evaluate bus angles
    data.bus_angles[data.valid_ix] =
        data.power_network_matrix.K \ power_injection[data.valid_ix]
    # evaluate flows
    my_mul_mt!(data.branch_flow_values, data.aux_network_matrix.data, data.bus_angles)
    return
end

# SINGLE PERIOD: method based on PTDF matrix
function _solve_powerflows_single!(
    data::PTDFPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    # evaluate flows
    my_mul_mt!(data.branch_flow_values, data.power_network_matrix.data, power_injection)
    # evaluate bus angles
    p_inj = power_injection[data.valid_ix]
    data.bus_angles[data.valid_ix] = data.aux_network_matrix.K \ p_inj
    return
end

# SINGLE PERIOD: method based on Virtual PTDF
function _solve_powerflows_single!(
    data::vPTDFPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    # evaluate flows (next line evaluates both both PTDF rows and line flows)
    my_mul_mt!(data.branch_flow_values, data.power_network_matrix, power_injection)
    # evaluate bus angles
    p_inj = power_injection[data.valid_ix]
    data.bus_angles[data.valid_ix] = data.aux_network_matrix.K \ p_inj
    return
end

# MULTI PERIOD: method based on ABA and BA matrices
function _solve_powerflows_multi!(
    data::ABAPowerFlowData,
)
    # get net injections
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    # save angles and power flows
    data.bus_angles[data.valid_ix, :] .=
        data.power_network_matrix.K \ @view power_injection[data.valid_ix, :]
    data.branch_flow_values .= data.aux_network_matrix.data' * data.bus_angles
    return
end

# MULTI PERIOD: method based on PTDF matrix
function _solve_powerflows_multi!(
    data::PTDFPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    # evaluate flows
    data.branch_flow_values .= data.power_network_matrix.data' * power_injection
    # evaluate bus angles
    p_inj = power_injection[data.valid_ix, :]
    data.bus_angles[data.valid_ix, :] .= data.aux_network_matrix.K \ p_inj
    return
end

# MULTI PERIOD: method based on Virtual PTDF
function _solve_powerflows_multi!(
    data::vPTDFPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    for i in axes(power_injection, 2)
        # evaluate flows (next line evaluates both both PTDF rows and line flows)
        data.branch_flow_values[:, i] .=
            my_mul_mt(data.power_network_matrix, power_injection[:, i])
        # evaluate bus angles
        p_inj = power_injection[data.valid_ix, i]
        data.bus_angles[data.valid_ix, i] .= data.aux_network_matrix.K \ p_inj
    end
    return
end

"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::PTDFPowerFlowData`:
        PTDFPowerFlowData structure containig all the information related to the system power flow
"""
function solve_powerflow!(
    data::PTDFPowerFlowData,
)
    if length(data.timestep_map) == 1
        _solve_powerflows_single!(data::PTDFPowerFlowData)
    else
        _solve_powerflows_multi!(data::PTDFPowerFlowData)
    end
    return
end

"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::vPTDFPowerFlowData`:
        vPTDFPowerFlowData structure containig all the information related to the system power flow
"""
function solve_powerflow!(
    data::vPTDFPowerFlowData,
)
    if length(data.timestep_map) == 1
        _solve_powerflows_single!(data::vPTDFPowerFlowData)
    else
        _solve_powerflows_multi!(data::vPTDFPowerFlowData)
    end
    return
end

# TODO: solve just for some lines with vPTDF

"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::ABAPowerFlowData`:
        ABAPowerFlowData structure containig all the information related to the system power flow
"""
# DC flow: ABA and BA case
function solve_powerflow!(
    data::ABAPowerFlowData,
)
    if length(data.timestep_map) == 1
        _solve_powerflows_single!(data::ABAPowerFlowData)
    else
        _solve_powerflows_multi!(data::ABAPowerFlowData)
    end
    return
end

function solve_powerflow(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(PTDFDCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys)
end

function solve_powerflow(
    ::DCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(DCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys)
end

function solve_powerflow(
    ::vPTDFDCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(vPTDFDCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys)
end
