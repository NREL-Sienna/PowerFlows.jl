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
    parallel::Bool,
)
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    matrix_data = data.power_network_matrix.data
    # Julia's Dense Matrix Multiplication is already multitheaded via OpenBlas
    LinearAlgebra.mul!(data.branch_flow_values, matrix_data, power_injection)
    valid_ix = data.valid_ix
    p_inj = view(power_injection, valid_ix)
    LinearAlgebra.ldiv!(data.bus_angle[valid_ix], data.aux_network_matrix.K, p_inj)
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
    parallel = false,
)
    data = PowerFlowData(DCPowerFlow(), sys)
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    matrix_data = data.power_network_matrix.K
    aux_network_matrix = data.aux_network_matrix

    # pending make this parallel and non-allocating
    data.aux_network_matrix.ref_bus_positions
    valid_ix = setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
    data.bus_angle[valid_ix] = matrix_data \ power_injection[valid_ix]
    if parallel
        my_mul_mt!(data.branch_flow_values, aux_network_matrix.data, data.bus_angle)
    else
        LinearAlgebra.mul!(
            data.branch_flow_values,
            aux_network_matrix.data,
            power_injection,
        )
    end
    return data
end
