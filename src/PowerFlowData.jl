struct PowerFlowData
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{Int, Int}
    bus_activepower_injection::Vector{Float64}
    bus_reactivepower_injection::Vector{Float64}
    bus_type::Vector{PSY.BusTypes}
    bus_magnitude::Vector{Float64}
    bus_angle::Vector{Float64}
    branch_flow_values::Vector{Float64}
    network_matrix::PNM.PowerNetworkMatrix
end

function PowerFlowData(::ACPowerFlow, sys::PSY.System)
    bus_ix = Dict{Int, Int}()
    branch_ix = Dict{Int, Int}()
    bus_activepower_injection = Vector{Float64}()
    bus_reactivepower_injection = Vector{Float64}()
    bus_type = Vector{PSY.BusTypes}()
    bus_magnitude = Vector{Float64}()
    bus_angle = Vector{Float64}()
    branch_activepower_flow_values = Vector{Float64}()
    branch_reactivepower_flow_values = Vector{Float64}()
    network_matrix = PNM.PowerNetworkMatrix(sys)
    for (ix, bus) in enumerate(PSY.get_components(PSY.Bus, sys))
        bus_ix[PSY.get_number(bus)] = ix
        bus_p_injection[ix] = PSY.get_active_power_injection(bus)
        bus_q_injection[ix] = PSY.get_reactive_power_injection(bus)
        bus_type[ix] = PSY.get_bustype(bus)
        bus_magnitude[ix] = PSY.get_magnitude(bus)
        bus_angle[ix] = PSY.get_angle(bus)
    end
    for (ix, branch) in enumerate(PSY.get_components(PSY.Branch, sys))
        branch_ix[PSY.get_name(branch)] = ix
        branch_flow_values[ix] = PSY.get_flow(branch)
    end
    return PowerFlowData(
        bus_ix,
        branch_ix,
        bus_p_injection,
        bus_q_injection,
        bus_type,
        bus_magnitude,
        bus_angle,
        branch_flow_values,
        network_matrix,
    )
end

function PowerFlowData(
    ::DCPowerFlow,
    sys::PSY.System,
    power_network_matrix::PNM.PowerNetworkMatrix,
)
    bus_ix = Dict{Int, Int}()
    branch_ix = Dict{Int, Int}()
    bus_p_injection = Vector{Float64}()
    bus_q_injection = Vector{Float64}()
    bus_type = Vector{PSY.BusTypes}()
    bus_magnitude = Vector{Float64}()
    bus_angle = Vector{Float64}()
    branch_flow_values = Vector{Float64}()
    network_matrix = power_network_matrix
    for (ix, bus) in enumerate(PSY.get_components(PSY.Bus, sys))
        bus_ix[PSY.get_number(bus)] = ix
        bus_p_injection[ix] = PSY.get_active_power_injection(bus)
        bus_type[ix] = PSY.get_bustype(bus)
        bus_angle[ix] = PSY.get_angle(bus)
    end
    for (ix, branch) in enumerate(PSY.get_components(PSY.Branch, sys))
        branch_ix[PSY.get_name(branch)] = ix
        branch_flow_values[ix] = PSY.get_active_power_flow(branch)
    end
    return PowerFlowData(
        bus_ix,
        branch_ix,
        bus_p_injection,
        bus_q_injection,
        bus_type,
        bus_magnitude,
        bus_angle,
        branch_flow_values,
        network_matrix,
    )
end
