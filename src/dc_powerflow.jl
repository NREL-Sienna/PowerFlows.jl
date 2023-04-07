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
    power_network_matrix::PNM.PTDF,
    aux_network_matrix::PNM.ABA_Matrix, data = PowerFlow

    ref_buses = aux_network_matrix.ref_bus_positions
    # ! setdiff still present since ptdf has also the columns related to the reference buses
    matrix_data = @view power_network_matrix.data[:, setdiff(1:end, ref_buses)]
    if parallel
        my_mul_mt!(data.branch_flow_values, matrix_data, power_injection)
    else
        my_mul_single!(data.branch_flow_values, matrix_data, power_injection)
    end

    ldiv!(data.bus_angle, aux_network_matrix.K, power_injection)
    return
end

function solve_powerflow!(
    ::DCPowerFlow,
    sys::PSY.System;
    parallel = false,
)
    data::PowerFlowData,
    power_network_matrix::PNM.ABA_Matrix,
    aux_network_matrix::PNM.BA_Matrix,
    power_injection::Vector{Float64}
    data.bus_angle[:] .= power_network_matrix.K \ power_injection
    data.branch_flow_values[:] .= aux_network_matrix.data * data.bus_angle
    return
end
