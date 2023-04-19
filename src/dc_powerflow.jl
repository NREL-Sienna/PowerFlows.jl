const PTDFPowerFlowData = PowerFlowData{
    PNM.PTDF{
        Tuple{Vector{String}, Vector{Int64}},
        Tuple{Dict{String, Int64}, Dict{Int64, Int64}},
        Matrix{Float64},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `::DCPowerFlow`:
        type of power flow analysis
- `data::PowerFlowData`:
        PowerFlowData structure containig all the information related to the system power flow
"""
function solve_powerflow!(
    data::PTDFPowerFlowData,
    parallel::Bool, # TODO still to implement
)
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    matrix_data = data.power_network_matrix.data
    if parallel
        my_mul_mt!(data.branch_flow_values, matrix_data, power_injection)
    else
        LinearAlgebra.mul!(data.branch_flow_values, matrix_data, power_injection)
    end

    #if parallel
    # to be added
    #else
    valid_ix =
        setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
    p_inj = power_injection[valid_ix]
    data.bus_angle[valid_ix] = data.aux_network_matrix.K \ p_inj
    #end
    return data
end

function solve_powerflow!(
    data::DCPowerFlow,
    parallel::Bool, # TODO still to implement
)

    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    matrix_data = data.power_network_matrix.K
    aux_network_matrix = data.aux_network_matrix

    # pending make this parallel and non-allocating
    valid_ix = setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
    data.bus_angle[valid_ix] = matrix_data \ power_injection[valid_ix]

    # evaluate
    # if parallel
    my_mul_mt_1!(data.branch_flow_values, aux_network_matrix.data, data.bus_angle)
    # else
    #     LinearAlgebra.mul!(
    #         data.branch_flow_values,
    #         aux_network_matrix.data,
    #         power_injection,
    #     )
    # end
    return data

end

function solve_powerflow!(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
    parallel = false,
)
    data = PowerFlowData(PTDFDCPowerFlow(), sys)
    return solve_powerflow!(data, parallel)
end

function solve_powerflow!(
    ::DCPowerFlow,
    sys::PSY.System;
    parallel = false,   # TODO still to implement
)
    data = PowerFlowData(DCPowerFlow(), sys)
    return solve_powerflow!(data, parallel)
end
