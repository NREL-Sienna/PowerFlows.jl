struct PowerFlowData{M, N <: PNM.PowerNetworkMatrix}
    n_buses::Int
    n_branches::Int
    ref_buses::Vector{Int}
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{Union{String, Int}, Int}
    bus_components::Dict{Int, Vector{PSY.StaticInjection}}
    bus_activepower_injection::Vector{Float64}
    bus_reactivepower_injection::Vector{Float64}
    bus_type::Vector{Any} # ! Vector{PSY.BusTypes}
    bus_magnitude::Vector{Float64}
    bus_angle::Vector{Float64}
    branch_flow_values::Vector{Float64}
    power_network_matrix::M
    aux_network_matrix::N
end

function PowerFlowData(::ACPowerFlow, sys::PSY.System)
    bus_lookup = Dict{Int, Int}()
    branch_lookup = Dict{String, Int}()
    bus_activepower_injection = Vector{Float64}()
    bus_reactivepower_injection = Vector{Float64}()
    bus_type = Vector{PSY.BusTypes}()
    bus_magnitude = Vector{Float64}()
    bus_angle = Vector{Float64}()
    branch_activepower_flow_values = Vector{Float64}()
    branch_reactivepower_flow_values = Vector{Float64}()
    network_matrix = PNM.PowerNetworkMatrix(sys)
    for (ix, bus) in enumerate(PSY.get_components(PSY.Bus, sys))
        bus_lookup[PSY.get_number(bus)] = ix
        bus_p_injection[ix] = PSY.get_active_power_injection(bus)
        bus_q_injection[ix] = PSY.get_reactive_power_injection(bus)
        bus_type[ix] = PSY.get_bustype(bus)
        bus_magnitude[ix] = PSY.get_magnitude(bus)
        bus_angle[ix] = PSY.get_angle(bus)
    end
    for (ix, branch) in enumerate(PSY.get_components(PSY.Branch, sys))
        branch_lookup[PSY.get_name(branch)] = ix
        branch_flow_values[ix] = PSY.get_flow(branch)
    end
    return PowerFlowData(
        bus_lookup,
        branch_lookup,
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

    # get number of slack buses
    if typeof(power_network_matrix) == PNM.ABA_Matrix
        ref_buses = power_network_matrix.ref_bus_positions
        n_ref_buses = length(power_network_matrix.ref_bus_positions)
    else
        ref_buses = aux_network_matrix.ref_bus_positions
        n_ref_buses = length(aux_network_matrix.ref_bus_positions)
    end

    # get number of buses and branches
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    n_branches = length(PSY.get_components(PSY.ACBranch, sys))  # ! PSY.Branch or PSY.ACBranch ???

    # Initizalize data
    # ! IMPORTANT: data related to ref buses is avoided
    bus_lookup = Dict{Int, Int}()
    branch_lookup = Dict{Union{String, Int}, Int}()
    bus_components = Dict{Int, Vector{PSY.StaticInjection}}()
    bus_activepower_injection = Vector{Float64}(undef, n_buses - n_ref_buses)
    bus_reactivepower_injection = Vector{Float64}(undef, n_buses - n_ref_buses)
    bus_type = Vector{Any}(undef, n_buses - n_ref_buses) # !!! should be PSY.BusTypes !!!
    bus_magnitude = Vector{Float64}(undef, n_buses - n_ref_buses)
    bus_angle = Vector{Float64}(undef, n_buses - n_ref_buses)
    branch_flow_values = Vector{Float64}(undef, n_branches)

    # get all the buse values (injection, angle, voltage, etc...)
    # ! IMPORTANT: no info related to ref buses
    for (ix, bus) in enumerate(PSY.get_components(PSY.Bus, sys))
        if ix ∉ ref_buses
            check_ = sum(ix .> ref_buses)
            bus_lookup[PSY.get_number(bus)] = ix
            bus_components[ix] = [d for d in PSY.get_components(PSY.StaticInjection, sys) if PSY.get_bus(d) == bus]
            bus_activepower_injection[ix - check_] = get_active_power_injection(bus_components[ix])
            bus_type[ix - check_] = PSY.get_bustype(bus)
            bus_angle[ix - check_] = PSY.get_angle(bus)
        end
    end

    # get lines' flows # ! PSY.Branch or PSY.ACBranch ???
    for (ix, branch) in enumerate(PSY.get_components(PSY.ACBranch, sys))
        branch_lookup[PSY.get_name(branch)] = ix
        # active power flow saved in branch
        branch_flow_values[ix] = PSY.get_active_power_flow(branch)
    end

    return PowerFlowData(
        n_buses,
        n_branches,
        ref_buses,
        bus_lookup,
        branch_lookup,
        bus_components,
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
