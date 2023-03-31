"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `::DCPowerFlow`:
        type of power flow analysis
- `data::PowerFlowData`:
        PowerFlowData structure containig all the information related to the system power flow
"""
function solve_powerflow(
    ::DCPowerFlow,
    data::PowerFlowData,
    bus_activepower_injection::Vector{Float64})
    data.branch_flow_values[:] .= _get_dc_flows(data.power_network_matrix,
        data.aux_power_network_matrix,
        power_injection)
    return data.branch_flow_values
end

function _get_dc_flows(
    power_network_matrix::PNM.PTDF,
    aux_power_network_matrix::Nothing,
    power_injection::Vector{Float64})
    data.branch_flow_values[:] .= power_network_matrix.data * power_injection
end

function _get_dc_flows(
    power_network_matrix::PNM.ABA_Matrix,
    aux_power_network_matrix::PNM.BA_Matrix,
    power_injection::Vector{Float64})
    thetas = power_network_matrix.K \ power_injection[setdiff(1:end, BA.ref_bus_positions)]
    data.branch_flow_values[:] .= aux_power_network_matrix.data * thetas
end

# get active power injected in a certain bus
function get_active_power_injection(sys::PSY.System, bus::PSY.Bus)
    # get the components connected to the considered bus
    comp_ =
        [d for d in PSY.get_components(PSY.StaticInjection, sys) if PSY.get_bus(d) == bus]
    # return the total injection at the cosndiered bus
    if length(comp_) > 0
        return sum(PSY.get_active_power(d) for d in comp_)
    else
        return 0
    end
end
