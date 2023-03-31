"""
Evaluates the power flowing on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `::DCPowerFlow`:
        type of power flow analysis
- `data::PowerFlowData`:
        PowerFlowData structure containig all the information related to the system power flow
"""
function solve_powerflow!(
    ::DCPowerFlow,
    data::PowerFlowData,
    power_injection::Vector{Float64})
    out_ = _get_dc_flows(data.network_matrix,
        data.aux_network_matrix,
        power_injection)
    data.bus_angle[:] .= out_[1]
    data.branch_flow_values[:] .= out_[2]
end

function _get_dc_flows(
    power_network_matrix::PNM.PTDF,
    aux_network_matrix::PNM.ABA_Matrix,
    power_injection::Vector{Float64})
    thetas = zeros(length(power_injection))
    thetas[setdiff(1:end, aux_network_matrix.ref_bus_positions)] =
        aux_network_matrix.K \
        power_injection[setdiff(1:end, aux_network_matrix.ref_bus_positions)]
    flows = power_network_matrix.data * power_injection
    out_ = (thetas, flows)

    return out_
end

function _get_dc_flows(
    power_network_matrix::PNM.ABA_Matrix,
    aux_network_matrix::PNM.BA_Matrix,
    power_injection::Vector{Float64})
    thetas = zeros(length(power_injection))
    thetas[setdiff(1:end, aux_network_matrix.ref_bus_positions)] =
        power_network_matrix.K \
        power_injection[setdiff(1:end, aux_network_matrix.ref_bus_positions)]
    flows =
        aux_network_matrix.data *
        thetas[setdiff(1:end, aux_network_matrix.ref_bus_positions)]
    out_ = (thetas, flows)

    return out_
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
