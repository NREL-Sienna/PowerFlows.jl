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

"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `::DCPowerFlow`:
        type of power flow analysis
- `data::PowerFlowData`:
        PowerFlowData structure containig all the information related to the system power flow
"""
# TODO consider adding argument "parallel::Bool", still to implement
# TODO missing iteration over columns
function solve_powerflow!(
    data::PTDFPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    matrix_data = data.power_network_matrix.data

    # evaluate flows
    # for vPTDFDCPowerFlow case, evaluates both PTDF rows and line flows
    my_mul_mt!(data.branch_flow_values, matrix_data, power_injection)
    p_inj = power_injection[data.valid_ix]
    # evaluate bus angles
    data.bus_angles[data.valid_ix] = data.aux_network_matrix.K \ p_inj

    return
end

# TODO consider adding argument "parallel::Bool", still to implement
# TODO missing iteration over columns
function solve_powerflow!(
    data::vPTDFPowerFlowData,
)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    matrix_data = data.power_network_matrix

    # evaluate flows
    # for vPTDFDCPowerFlow case, evaluates both PTDF rows and line flows
    my_mul_mt!(data.branch_flow_values, matrix_data, power_injection)
    p_inj = power_injection[data.valid_ix]
    # evaluate bus angles
    data.bus_angles[data.valid_ix] = data.aux_network_matrix.K \ p_inj

    return
end

# TODO consider adding argument "parallel::Bool", still to implement
function solve_powerflow_for!(
    data::ABAPowerFlowData,
)
    # get matrices
    matrix_data = data.power_network_matrix.K
    aux_network_matrix = data.aux_network_matrix
    for i in axes(data.bus_activepower_injection, 2)
        # get net power injections
        power_injection = data.bus_activepower_injection[:, i] - data.bus_activepower_withdrawals[:, i]
        # evaluate bus angles
        data.bus_angles[data.valid_ix, i] .= matrix_data \ power_injection[data.valid_ix]
        # evaluate flows
        data.branch_flow_values[:, i] .= my_mul_mt(aux_network_matrix.data, @view data.bus_angles[:, i])
    end
    return
end

# ! reference function for ABAPowerFlowData type
function solve_powerflow!(
    data::ABAPowerFlowData,
)
    # get matrices
    matrix_data = data.power_network_matrix.K
    aux_network_matrix = data.aux_network_matrix
    # get net power injections
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    # evaluate bus angles
    data.bus_angles[data.valid_ix] = matrix_data \ power_injection[data.valid_ix]
    # evaluate flows
    my_mul_mt!(data.branch_flow_values, aux_network_matrix.data, data.bus_angles)
    return
end

# TODO consider adding argument "parallel::Bool", still to implement
function solve_powerflow(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(PTDFDCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys) # TODO missing iteration over columns
end

# TODO consider adding argument "parallel::Bool", still to implement
function solve_powerflow(
    ::DCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(DCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys) # TODO missing iteration over columns
end

# TODO consider adding argument "parallel::Bool", still to implement
function solve_powerflow(
    ::vPTDFDCPowerFlow,
    sys::PSY.System;
)
    data = PowerFlowData(vPTDFDCPowerFlow(), sys)
    solve_powerflow!(data)
    return write_results(data, sys) # TODO missing iteration over columns
end
