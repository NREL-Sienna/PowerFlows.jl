
# TODO a little clunky: all of the form length(get_{arc/bus}_axis( )).
arc_count(::ACPowerFlow,
    power_network_matrix::PNM.PowerNetworkMatrix,
    ::Union{PNM.PowerNetworkMatrix, Nothing}) = length(PNM.get_arc_axis(power_network_matrix.branch_admittance_from_to))
bus_count(::ACPowerFlow,
    power_network_matrix::PNM.PowerNetworkMatrix,
    ::Union{PNM.PowerNetworkMatrix, Nothing}) = length(PNM.get_bus_axis(power_network_matrix))

arc_count(::Union{PTDFDCPowerFlow, vPTDFDCPowerFlow},
    power_network_matrix::PNM.PowerNetworkMatrix,
    ::Union{PNM.PowerNetworkMatrix, Nothing}) =
    length(PNM.get_arc_axis(power_network_matrix))
bus_count(::Union{PTDFDCPowerFlow, vPTDFDCPowerFlow},
    power_network_matrix::PNM.PowerNetworkMatrix,
    ::Union{PNM.PowerNetworkMatrix, Nothing}) =
    length(PNM.get_bus_axis(power_network_matrix))

arc_count(::DCPowerFlow,
    power_network_matrix::PNM.PowerNetworkMatrix,
    aux_network_matrix::Union{PNM.PowerNetworkMatrix, Nothing}) =
    length(PNM.get_arc_axis(aux_network_matrix))
bus_count(::DCPowerFlow,
    power_network_matrix::PNM.PowerNetworkMatrix,
    aux_network_matrix::Union{PNM.PowerNetworkMatrix, Nothing}) =
    length(PNM.get_bus_axis(aux_network_matrix))

function make_powerflow_data(
    pf::T,
    power_network_matrix::M,
    aux_network_matrix::N,
    n_timesteps::Int;
    timestep_names::Vector{String} = String[],
    neighbors = Vector{Set{Int}}(),
) where {
    T <: PowerFlowEvaluationModel,
    M <: PNM.PowerNetworkMatrix,
    N <: Union{PNM.PowerNetworkMatrix, Nothing},
}
    if n_timesteps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:n_timesteps]
        elseif length(timestep_names) != n_timesteps
            error("timestep_names field must have same length as n_timesteps")
        end
    end
    timestep_map = Dict(zip([i for i in 1:n_timesteps], timestep_names))

    n_buses = bus_count(pf, power_network_matrix, aux_network_matrix)
    n_arcs = arc_count(pf, power_network_matrix, aux_network_matrix)
    calculate_loss_factors = get_calculate_loss_factors(pf)
    calculate_voltage_stability_factors = get_calculate_voltage_stability_factors(pf)
    if !isnothing(get_slack_participation_factors(pf))
        empty_slack_participation_factors = Dict{Tuple{DataType, String}, Float64}[]
    else
        empty_slack_participation_factors = nothing
    end
    return PowerFlowData(
        zeros(n_buses, n_timesteps), # bus_activepower_injection
        zeros(n_buses, n_timesteps), # bus_reactivepower_injection
        zeros(n_buses, n_timesteps), # bus_activepower_withdrawals
        zeros(n_buses, n_timesteps), # bus_reactivepower_withdrawals
        zeros(n_buses, n_timesteps), # bus_activepower_constant_current_withdrawals
        zeros(n_buses, n_timesteps), # bus_reactivepower_constant_current_withdrawals
        zeros(n_buses, n_timesteps), # bus_activepower_constant_impedance_withdrawals
        zeros(n_buses, n_timesteps), # bus_reactivepower_constant_impedance_withdrawals
        fill((-Inf, Inf), (n_buses, n_timesteps)), # bus_reactivepower_bounds
        empty_slack_participation_factors, # generator_slack_participation_factors
        spzeros(n_buses, n_timesteps), # bus_slack_participation_factors
        fill(PSY.ACBusTypes.PQ, (n_buses, n_timesteps)), # bus_type
        ones(n_buses, n_timesteps), # bus_magnitude
        zeros(n_buses, n_timesteps), # bus_angles
        zeros(n_arcs, n_timesteps), # arc_activepower_flow_from_to
        zeros(n_arcs, n_timesteps), # arc_reactivepower_flow_from_to
        zeros(n_arcs, n_timesteps), # arc_activepower_flow_to_from
        zeros(n_arcs, n_timesteps), # arc_reactivepower_flow_to_from
        timestep_map,
        power_network_matrix,
        aux_network_matrix,
        neighbors,
        falses(n_timesteps), # converged
        calculate_loss_factors ? zeros(n_buses, n_timesteps) : nothing, # loss_factors
        calculate_loss_factors,
        calculate_voltage_stability_factors ? zeros(n_buses, n_timesteps) : nothing, # voltage_stability_factors
        calculate_voltage_stability_factors,
    )
end
