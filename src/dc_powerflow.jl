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
        data.branch_flow_values[:] = matrix_data*power_injection
    end

    # pending make this parallel
    # Not working yet BA doesn't have the same dimensions as PowerInjection
    #data.bus_angle[:] = data.aux_network_matrix.data \ power_injection
    return
end

function solve_powerflow!(
    ::DCPowerFlow,
    sys::PSY.System;
    parallel = false,
)
    data = PowerFlowData(DCPowerFlow(), sys)

    # pending make this parallel
    data.bus_angle[:] .= power_network_matrix.K \ power_injection
    if parallel
        my_mul_mt!(data.branch_flow_values, aux_network_matrix.data, data.bus_angle)
    else
        data.branch_flow_values[:] = aux_network_matrix.data * data.bus_angle
    end
    return
end
