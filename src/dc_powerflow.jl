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
    if typeof(data.power_network_matrix) <: PNM.PTDF &&
       typeof(data.aux_network_matrix) <: PNM.ABA_Matrix
        data.bus_angle[:] .= data.aux_network_matrix.K \ power_injection
        ref_buses = data.aux_network_matrix.ref_bus_positions
        # ! setdiff still present since ptdf has also the columns related to the reference buses
        data.branch_flow_values[:] .=
            data.power_network_matrix.data[:, setdiff(1:end, ref_buses)] * power_injection
    elseif typeof(data.power_network_matrix) <: PNM.ABA_Matrix &&
           typeof(data.aux_network_matrix) <: PNM.BA_Matrix
        data.bus_angle[:] .= data.power_network_matrix.K \ power_injection
        data.branch_flow_values[:] .= data.aux_network_matrix.data * data.bus_angle[:]
    end
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
