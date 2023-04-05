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
    power_network_matrix::PNM.PTDF,
    aux_network_matrix::PNM.ABA_Matrix,
    power_injection::Vector{Float64}
    )

    data.bus_angle[:] .= aux_network_matrix.K \ power_injection
    ref_buses = aux_network_matrix.ref_bus_positions
    # ! setdiff still present since ptdf has also the columns related to the reference buses
    matrix_data = @view power_network_matrix.data[:, setdiff(1:end, ref_buses)]
    data.branch_flow_values[:] .= matrix_data * power_injection
end

function solve_powerflow!(
    ::DCPowerFlow,
    data::PowerFlowData,
    power_network_matrix::PNM.ABA_Matrix,
    aux_network_matrix::PNM.BA_Matrix,
    power_injection::Vector{Float64}
    )

    data.bus_angle[:] .= power_network_matrix.K \ power_injection
    data.branch_flow_values[:] .= aux_network_matrix.data * data.bus_angle
end

# get active power injected in a certain bus
function get_active_power_injection(bus_components::Vector{T}) where T <: PSY.StaticInjection
    # get the components connected to the considered bus
    pos_comp = [d for d in filter(x -> !isa(x, PSY.ElectricLoad), bus_components)]
    neg_comp = [d for d in filter(x -> isa(x, PSY.ElectricLoad), bus_components)]
    # return the total injection at the cosndiered bus
    return _sum_active_power(pos_comp) - _sum_active_power(neg_comp)
end

function _sum_active_power(bus_components::Vector{T}) where T <: PSY.StaticInjection
    vec_ = [PSY.get_active_power(d) for d in bus_components if typeof(d) != PSY.FixedAdmittance]
    if length(vec_) > 0
        return sum(vec_)
    else
        return 0
    end
end