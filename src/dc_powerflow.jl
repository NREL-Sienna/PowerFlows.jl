"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `::DCPowerFlow`:
        type of power flow analysis
- `data::PowerFlowData`:
        PowerFlowData structure containig all the information related to the system power flow
"""
function solve_powerflow!(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
    parallel = false,
)
    data = PowerFlowData(PTDFDCPowerFlow(), sys)
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    matrix_data = data.power_network_matrix.data
    if parallel
        my_mul_mt!(data.branch_flow_values, matrix_data, power_injection)
    else
        LinearAlgebra.mul!(data.branch_flow_values, matrix_data, power_injection)
    end

    if parallel
        # to be added
    else
        valid_ix =
            setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
        p_inj = power_injection[valid_ix]
        data.bus_angle[valid_ix] = data.aux_network_matrix.K \ p_inj
    end
    return data
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
