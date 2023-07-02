# ! NOTE: bus_magnitude saved as a matrix like for angles and line flows. Is it correct? 
# ! or just keep a vector?

struct PowerFlowData{M <: PNM.PowerNetworkMatrix, N, S <: Union{String, Char}}
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{String, Int}
    bus_activepower_injection::Union{Vector{Float64}, Matrix{Float64}}
    bus_reactivepower_injection::Union{Vector{Float64}, Matrix{Float64}}
    bus_activepower_withdrawals::Union{Vector{Float64}, Matrix{Float64}}
    bus_reactivepower_withdrawals::Union{Vector{Float64}, Matrix{Float64}}
    bus_type::Vector{PSY.BusTypes}
    bus_magnitude::Union{Vector{Float64}, Matrix{Float64}}
    bus_angles::Union{Vector{Float64}, Matrix{Float64}}
    branch_flow_values::Union{Vector{Float64}, Matrix{Float64}}
    timestep_map::Dict{Int, S}
    valid_ix::Vector{Int}
    power_network_matrix::M
    aux_network_matrix::N
end

# SINGLE PERIOD: AC Power Flow Data
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
        zeros(n_branches, 1),
        Dict(zip([1], "1")),
        setdiff(1:n_buses, ref_bus_positions),
        power_network_matrix,
        nothing,
    )
end

# SINGLE PERIOD: DC Power Flow Data based on ABA and BA matrices
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
    bus_angles = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
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
        bus_angles,
        zeros(n_branches),
        Dict(zip([1], "1")),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )

end

# SINGLE PERIOD: DC Power Flow Data based on PTDF matrix
function PowerFlowData(::PTDFDCPowerFlow, sys::PSY.System)

    # get the network matrices
    power_network_matrix = PNM.PTDF(sys)
    aux_network_matrix = PNM.ABA_Matrix(sys; factorize = true)

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 1))
    n_branches = length(axes(power_network_matrix, 2))

    bus_lookup = power_network_matrix.lookup[1]
    branch_lookup = power_network_matrix.lookup[2]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
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
        bus_angles,
        zeros(Float64, n_branches),
        Dict(zip([1], "1")),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end

# SINGLE PERIOD: DC Power Flow Data based on virutual PTDF matrix
function PowerFlowData(::vPTDFDCPowerFlow, sys::PSY.System)

    # get the network matrices
    power_network_matrix = PNM.VirtualPTDF(sys) # evaluates an empty virtual PTDF
    aux_network_matrix = PNM.ABA_Matrix(sys; factorize = true)

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 2))
    n_branches = length(axes(power_network_matrix, 1))

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup = power_network_matrix.lookup[1]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
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
        bus_angles,
        zeros(Float64, n_branches),
        Dict(zip([1], "1")),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end

# TODO -> MULTI PERIOD: AC Power Flow Data

# MULTI PERIOD: DC Power Flow Data based on ABA and BA matrices
function PowerFlowData(
    ::DCPowerFlow,
    sys::PSY.System,
    timesteps::Int,
    timestep_names::Vector{String} = String[])

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if timesteps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:timesteps]
        elseif  length(timestep_names) != timesteps
            error("timestep_names field must have same length as timesteps")
        end
    end

    # get the network matrices
    power_network_matrix = PNM.ABA_Matrix(sys; factorize = true)
    aux_network_matrix = PNM.BA_Matrix(sys)

    # get number of buses and branches
    n_buses = length(axes(aux_network_matrix, 2))
    n_branches = length(axes(aux_network_matrix, 1))

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup = aux_network_matrix.lookup[1]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
    end

    # define injection vectors related to the first timestep
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

    # initialize data
    init_1 = zeros(n_buses, timesteps)
    init_2 = zeros(n_branches, timesteps)

    # define fields as matrices whose number of columns is eqault to the number of timesteps
    bus_activepower_injection_1 = deepcopy(init_1)
    bus_reactivepower_injection_1 = deepcopy(init_1)
    bus_activepower_withdrawals_1 = deepcopy(init_1)
    bus_reactivepower_withdrawals_1 = deepcopy(init_1)
    bus_magnitude_1 = deepcopy(init_1)
    bus_angles_1 = deepcopy(init_1)
    branch_flow_values_1 = deepcopy(init_2)

    # initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= ones(Float64, n_buses)
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        bus_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:timesteps], timestep_names)),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end

# MULTI PERIOD: DC Power Flow Data with PTDF matrix
function PowerFlowData(
    ::PTDFDCPowerFlow,
    sys::PSY.System,
    timesteps::Int,
    timestep_names::Vector{String} = String[])

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if timesteps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:timesteps]
        elseif  length(timestep_names) != timesteps
            error("timestep_names field must have same length as timesteps")
        end
    end

    # get the network matrices
    power_network_matrix = PNM.PTDF(sys)
    aux_network_matrix = PNM.ABA_Matrix(sys; factorize = true)

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 1))
    n_branches = length(axes(power_network_matrix, 2))

    bus_lookup = power_network_matrix.lookup[1]
    branch_lookup = power_network_matrix.lookup[2]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
    end

    # define injection vectors related to the first timestep
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

    # initialize data
    init_1 = zeros(n_buses, timesteps)
    init_2 = zeros(n_branches, timesteps)

    # define fields as matrices whose number of columns is eqault to the number of timesteps
    bus_activepower_injection_1 = deepcopy(init_1)
    bus_reactivepower_injection_1 = deepcopy(init_1)
    bus_activepower_withdrawals_1 = deepcopy(init_1)
    bus_reactivepower_withdrawals_1 = deepcopy(init_1)
    bus_magnitude_1 = deepcopy(init_1)
    bus_angles_1 = deepcopy(init_1)
    branch_flow_values_1 = deepcopy(init_2)

    # initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= ones(Float64, n_buses)
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        bus_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:timesteps], timestep_names)),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end

# MULTI PERIOD: DC Power Flow Data with virtual PTDF matrix
function PowerFlowData(
    ::vPTDFDCPowerFlow,
    sys::PSY.System,
    timesteps::Int,
    timestep_names::Vector{String} = String[])

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if timesteps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:timesteps]
        elseif  length(timestep_names) != timesteps
            error("timestep_names field must have same length as timesteps")
        end
    end

    # get the network matrices
    power_network_matrix = PNM.VirtualPTDF(sys) # evaluates an empty virtual PTDF
    aux_network_matrix = PNM.ABA_Matrix(sys; factorize = true)

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 2))
    n_branches = length(axes(power_network_matrix, 1))

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup = power_network_matrix.lookup[1]
    bus_type = Vector{PSY.BusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.BusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
    end

    # define injection vectors related to the first timestep
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

    # initialize data
    init_1 = zeros(n_buses, timesteps)
    init_2 = zeros(n_branches, timesteps)

    # define fields as matrices whose number of columns is eqault to the number of timesteps
    bus_activepower_injection_1 = deepcopy(init_1)
    bus_reactivepower_injection_1 = deepcopy(init_1)
    bus_activepower_withdrawals_1 = deepcopy(init_1)
    bus_reactivepower_withdrawals_1 = deepcopy(init_1)
    bus_magnitude_1 = deepcopy(init_1)
    bus_angles_1 = deepcopy(init_1)
    branch_flow_values_1 = deepcopy(init_2)

    # initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= ones(Float64, n_buses)
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        bus_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:timesteps], timestep_names)),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
    )
end