struct PowerFlowData{M, N <: PNM.PowerNetworkMatrix}
    n_buses::Int
    n_branches::Int
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{Union{String, Int}, Int}
    bus_activepower_injection::Vector{Float64}
    bus_reactivepower_injection::Vector{Float64}
    bus_type::Vector{Any} # ! Vector{PSY.BusTypes}
    bus_magnitude::Vector{Float64}
    bus_angle::Vector{Float64}
    branch_flow_values::Vector{Float64}
    network_matrix::M
    aux_network_matrix::N
end

function PowerFlowData(::ACPowerFlow, sys::PSY.System)
    bus_ix = Dict{Int, Int}()
    branch_ix = Dict{String, Int}()
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

# version with full PTDF
function PowerFlowData(
    ::DCPowerFlow,
    sys::PSY.System,
    power_network_matrix::Union{PNM.PTDF, PNM.ABA_Matrix},
    aux_network_matrix::Union{PNM.BA_Matrix, PNM.ABA_Matrix},
)

    # get number of buses and branches
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    n_branches = length(PSY.get_components(PSY.ACBranch, sys))  # ! PSY.Branch or PSY.ACBranch ???

    # Initizalize data
    bus_ix = Dict{Int, Int}()
    branch_ix = Dict{Union{String, Int}, Int}()
    bus_activepower_injection = Vector{Float64}(undef, n_buses)
    bus_reactivepower_injection = Vector{Float64}(undef, n_buses)
    bus_type = Vector{Any}(undef, n_buses) # !!! should be PSY.BusTypes !!!
    bus_magnitude = Vector{Float64}(undef, n_buses)
    bus_angle = Vector{Float64}(undef, n_buses)
    branch_flow_values = Vector{Float64}(undef, n_branches)

    # get all the buse values (injection, angle, voltage, etc...)
    for (ix, bus) in enumerate(PSY.get_components(PSY.Bus, sys))
        bus_ix[PSY.get_number(bus)] = ix
        bus_activepower_injection[ix] = get_active_power_injection(sys, bus)
        bus_type[ix] = PSY.get_bustype(bus)
        bus_angle[ix] = PSY.get_angle(bus)
    end

    # get lines' flows # ! PSY.Branch or PSY.ACBranch ???
    for (ix, branch) in enumerate(PSY.get_components(PSY.ACBranch, sys))
        branch_ix[PSY.get_name(branch)] = ix
        # active power flow saved in branch
        branch_flow_values[ix] = PSY.get_active_power_flow(branch)
    end

    return PowerFlowData(
        n_buses,
        n_branches,
        bus_ix,
        branch_ix,
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_type,
        bus_magnitude,
        bus_angle,
        branch_flow_values,
        power_network_matrix,
        aux_network_matrix,
    )
end
