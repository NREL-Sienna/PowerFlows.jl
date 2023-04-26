struct PowerFlowData{M <: PNM.PowerNetworkMatrix, N}
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{String, Int}
    bus_activepower_injection::Vector{Float64}
    bus_reactivepower_injection::Vector{Float64}
    bus_activepower_withdrawals::Vector{Float64}
    bus_reactivepower_withdrawals::Vector{Float64}
    bus_type::Vector{PSY.BusTypes}
    bus_magnitude::Vector{Float64}
    bus_angles::Vector{Float64}
    branch_flow_values::Vector{Float64}
    valid_ix::Vector{Int}
    power_network_matrix::M
    aux_network_matrix::N
end

function PowerFlowData(::ACPowerFlow, sys::PSY.System)
    power_network_matrix = PNM.Ybus(sys)

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 1))
    buses = PNM.get_buses(sys)
    ref_bus_positions = PNM.find_slack_positions(buses)

    branches = PNM.get_ac_branches(sys)
    n_branches = length(branches)

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup =
        Dict{String, Int}(PSY.get_name(b) => ix for (ix, b) in enumerate(branches))
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (ix, bus_no) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        bus_angles[ix] = PSY.get_angle(bus)
        bus_magnitude[ix] = PSY.get_magnitude(bus)
    end

    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    get_injections!(bus_activepower_injection, bus_reactivepower_injection, bus_lookup, sys)

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    get_withdrawals!(
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_lookup,
        sys,
    )

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_type,
        bus_magnitude,
        bus_angles,
        zeros(n_branches),
        setdiff(1:n_buses, ref_bus_positions),
        power_network_matrix,
        nothing,
    )
end

# version with full PTDF
function PowerFlowData(::DCPowerFlow, sys::PSY.System)
    power_network_matrix = PNM.ABA_Matrix(sys; factorize = true)
    aux_network_matrix = PNM.BA_Matrix(sys)
    # check the maps betwen the 2 matrices match

    # get number of buses and branches
    n_buses = length(axes(aux_network_matrix, 2))
    n_branches = length(axes(aux_network_matrix, 1))

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup = aux_network_matrix.lookup[1]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angle = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angle[ix] = 0.0
        else
            bus_angle[ix] = PSY.get_angle(bus)
        end
    end

    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    get_injections!(bus_activepower_injection, bus_reactivepower_injection, bus_lookup, sys)

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    get_withdrawals!(
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_lookup,
        sys,
    )

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_type,
        ones(Float64, n_buses),
        bus_angle,
        zeros(n_branches),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end

function PowerFlowData(::PTDFDCPowerFlow, sys::PSY.System)
    power_network_matrix = PNM.PTDF(sys)
    aux_network_matrix = PNM.ABA_Matrix(sys; factorize = true)
    # check the maps betwen the 2 matrices match

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 1))
    n_branches = length(axes(power_network_matrix, 2))

    bus_lookup = power_network_matrix.lookup[1]
    branch_lookup = power_network_matrix.lookup[2]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angle = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angle[ix] = 0.0
        else
            bus_angle[ix] = PSY.get_angle(bus)
        end
    end

    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    get_injections!(bus_activepower_injection, bus_reactivepower_injection, bus_lookup, sys)

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    get_withdrawals!(
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_lookup,
        sys,
    )

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_type,
        ones(Float64, n_buses),
        bus_angle,
        zeros(Float64, n_branches),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end

function PowerFlowData(::vPTDFDCPowerFlow, sys::PSY.System)
    power_network_matrix = PNM.VirtualPTDF(sys) # evaluates an empty virtual PTDF
    aux_network_matrix = PNM.ABA_Matrix(sys; factorize = true)

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 2))
    n_branches = length(axes(power_network_matrix, 1))

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup = power_network_matrix.lookup[1]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angle = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angle[ix] = 0.0
        else
            bus_angle[ix] = PSY.get_angle(bus)
        end
    end

    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    get_injections!(bus_activepower_injection, bus_reactivepower_injection, bus_lookup, sys)

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    get_withdrawals!(
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_lookup,
        sys,
    )

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_activepower_withdrawals,
        bus_reactivepower_withdrawals,
        bus_type,
        ones(Float64, n_buses),
        bus_angle,
        zeros(Float64, n_branches),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end
