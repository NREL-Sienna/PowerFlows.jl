"""
Structure containing all the data required for the evaluation of the power
flows and angles, as well as these ones.

# Arguments:
- `bus_lookup::Dict{Int, Int}`:
        dictionary linking the system's bus number with the rows of either
        "power_network_matrix" or "aux_network_matrix".
- `branch_lookup::Dict{String, Int}`:
        dictionary linking the branch name with the column name of either the
        "power_network_matrix" or "aux_network_matrix".
- `bus_activepower_injection::Matrix{Float64}`:
        "(b, t)" matrix containing the bus active power injection. b: number of
        buses, t: number of time period.
- `bus_reactivepower_injection::Matrix{Float64}`:
        "(b, t)" matrix containing the bus reactive power injection. b: number
        of buses, t: number of time period.
- `bus_activepower_withdrawals::Matrix{Float64}`:
        "(b, t)" matrix containing the bus reactive power withdrawals. b:
        number of buses, t: number of time period.
- `bus_reactivepower_withdrawals::Matrix{Float64}`:
        "(b, t)" matrix containing the bus reactive power withdrawals. b:
        number of buses, t: number of time period.
- `bus_reactivepower_bounds::Vector{Float64}`: Upper and Lower bounds for the reactive supply
        at each bus.
- `bus_type::Vector{PSY.ACBusTypes}`:
        vector containing type of buses present in the system, ordered
        according to "bus_lookup".
- `bus_magnitude::Matrix{Float64}`:
        "(b, t)" matrix containing the bus magnitudes, ordered according to
        "bus_lookup". b: number of buses, t: number of time period.
- `bus_angles::Matrix{Float64}`:
        "(b, t)" matrix containing the bus angles, ordered according to
        "bus_lookup". b: number of buses, t: number of time period.
- `branch_flow_values::Matrix{Float64}`:
        "(br, t)" matrix containing the power flows, ordered according to
        "branch_lookup". br: number of branches, t: number of time period.
- `timestep_map::Dict{Int, S}`:
        dictonary mapping the number of the time periods (corresponding to the
        column number of the previosly mentioned matrices) and their actual
        names.
- `valid_ix::Vector{Int}`:
        vector containing the indeces related to those buses that are not slack
        ones.
- `power_network_matrix::M`:
        matrix used for the evaluation of either the power flows or bus angles,
        depending on the method considered.
- `aux_network_matrix::N`:
        matrix used for the evaluation of either the power flows or bus angles,
        depending on the method considered.
- `neighbors::Vector{Set{Int}}`: Vector with the sets of adjacent buses.
- `x0::Matrix{Float64}`: Initial conditions for the AC PowerFlow.
"""
struct PowerFlowData{M <: PNM.PowerNetworkMatrix, N, S <: Union{String, String}}
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{String, Int}
    bus_activepower_injection::Matrix{Float64}
    bus_reactivepower_injection::Matrix{Float64}
    bus_activepower_withdrawals::Matrix{Float64}
    bus_reactivepower_withdrawals::Matrix{Float64}
    bus_reactivepower_bounds::Vector{Float64}
    bus_type::Vector{PSY.ACBusTypes}
    bus_magnitude::Matrix{Float64}
    bus_angles::Matrix{Float64}
    branch_flow_values::Matrix{Float64}
    timestep_map::Dict{Int, S}
    valid_ix::Vector{Int}
    power_network_matrix::M
    aux_network_matrix::N
    neighbors::Vector{Set{Int}}
end

# AC Power Flow Data
# TODO -> MULTI PERIOD: AC Power Flow Data
function _calculate_neighbors(
    Yb::PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
)
    I, J, V = SparseArrays.findnz(Yb.data)
    neighbors = [Set{Int}([i]) for i in 1:length(Yb.axes[1])]
    for nz in eachindex(V)
        push!(neighbors[I[nz]], J[nz])
        push!(neighbors[J[nz]], I[nz])
    end
    return neighbors
end

"""
Function for the definition of the PowerFlowData strucure given the System
data, number of time periods to consider and their names.
Calling this function will not evaluate the power flows and angles.
NOTE: use it for AC power flow computations.

# Arguments:
- `::ACPowerFlow`:
        use ACPowerFlow() to evaluate the AC OPF.
- `sys::PSY.System`:
        container storing the system data to consider in the PowerFlowData
        structure.
- `timesteps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "timesteps". Default
        value = String[].
- `check_connectivity::Bool`:
        Perform connectivity check on the network matrix. Default value = true.

WARNING: functions for the evaluation of the multi-period AC OPF still to be implemented.
"""
function PowerFlowData(
    ::ACPowerFlow,
    sys::PSY.System;
    timesteps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if timesteps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:timesteps]
        elseif length(timestep_names) != timesteps
            error("timestep_names field must have same length as timesteps")
        end
    end

    # get data for calculations
    power_network_matrix = PNM.Ybus(sys; check_connectivity = check_connectivity)

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 1))
    buses = PNM.get_buses(sys)
    ref_bus_positions = PNM.find_slack_positions(buses)

    branches = PNM.get_ac_branches(sys)
    n_branches = length(branches)

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup =
        Dict{String, Int}(PSY.get_name(b) => ix for (ix, b) in enumerate(branches))

    # TODO: bus_type might need to also be a Matrix since the type can change for a particular scenario
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.ACBusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
        bus_magnitude[ix] = PSY.get_magnitude(bus)
    end

    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    _get_injections!(
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_lookup,
        sys,
    )

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    _get_withdrawals!(
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
    bus_magnitude_1[:, 1] .= bus_magnitude
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        Vector{Float64}(),
        bus_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(1 => "1"),
        setdiff(1:n_buses, ref_bus_positions),
        power_network_matrix,
        nothing,
        _calculate_neighbors(power_network_matrix),
    )
end

# DC Power Flow Data based on ABA and BA matrices
"""
Function for the definition of the PowerFlowData strucure given the System
data, number of time periods to consider and their names.
Calling this function will not evaluate the power flows and angles.
NOTE: use it for DC power flow computations.

# Arguments:
- `::DCPowerFlow`:
        use DCPowerFlow() to store the ABA matrix as power_network_matrix and
        the BA matrix as aux_network_matrix.
- `sys::PSY.System`:
        container storing the system data to consider in the PowerFlowData
        structure.
- `timesteps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "timesteps". Default
        value = String[].
- `check_connectivity::Bool`:
        Perform connectivity check on the network matrix. Default value = true.
"""
function PowerFlowData(
    ::DCPowerFlow,
    sys::PSY.System;
    timesteps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if length(timestep_names) == 0
        timestep_names = [string(i) for i in 1:timesteps]
    elseif length(timestep_names) != timesteps
        error("timestep_names field must have same length as timesteps")
    end

    # get the network matrices
    power_network_matrix = PNM.ABA_Matrix(sys; factorize = true)
    aux_network_matrix = PNM.BA_Matrix(sys)

    # get number of buses and branches
    n_buses = length(axes(aux_network_matrix, 1))
    n_branches = length(axes(aux_network_matrix, 2))

    bus_lookup = aux_network_matrix.lookup[1]
    branch_lookup = aux_network_matrix.lookup[2]
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.ACBusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
        bus_magnitude[ix] = PSY.get_magnitude(bus)
    end

    # define injection vectors related to the first timestep
    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    _get_injections!(
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_lookup,
        sys,
    )

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    _get_withdrawals!(
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
    bus_magnitude_1 = zeros(n_buses, 1)
    bus_angles_1 = deepcopy(init_1)
    branch_flow_values_1 = deepcopy(init_2)

    # initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= bus_magnitude  # for DC case same value accross all timesteps
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        Vector{Float64}(),
        bus_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:timesteps], timestep_names)),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
        Vector{Set{Int}}(),
    )
end

# DC Power Flow Data with PTDF matrix
"""
Function for the definition of the PowerFlowData strucure given the System
data, number of time periods to consider and their names.
Calling this function will not evaluate the power flows and angles.
NOTE: use it for DC power flow computations.

# Arguments:
- `::PTDFDCPowerFlow`:
        use PTDFDCPowerFlow() to store the PTDF matrix as power_network_matrix
        and the ABA matrix as aux_network_matrix.
- `sys::PSY.System`:
        container storing the system data to consider in the PowerFlowData
        structure.
- `timesteps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "timesteps". Default
        value = String[].
"""
function PowerFlowData(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
    timesteps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if timesteps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:timesteps]
        elseif length(timestep_names) != timesteps
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
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.ACBusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
        bus_magnitude[ix] = PSY.get_magnitude(bus)
    end

    # define injection vectors related to the first timestep
    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    _get_injections!(
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_lookup,
        sys,
    )

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    _get_withdrawals!(
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
    bus_magnitude_1 = zeros(n_buses, 1)
    bus_angles_1 = deepcopy(init_1)
    branch_flow_values_1 = deepcopy(init_2)

    # initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= bus_magnitude  # for DC case same value accross all timesteps
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        Vector{Float64}(),
        bus_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:timesteps], timestep_names)),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
        Vector{Set{Int}}(),
    )
end

# DC Power Flow Data with virtual PTDF matrix
"""
Function for the definition of the PowerFlowData strucure given the System
data, number of time periods to consider and their names.
Calling this function will not evaluate the power flows and angles.
NOTE: use it for DC power flow computations.

# Arguments:
- `::PTDFDCPowerFlow`:
        use vPTDFDCPowerFlow() to store the Virtual PTDF matrix as
        power_network_matrix and the ABA matrix as aux_network_matrix.
- `sys::PSY.System`:
        container storing the system data to consider in the PowerFlowData
        structure.
- `timesteps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "timesteps". Default
        value = String[].
"""
function PowerFlowData(
    ::vPTDFDCPowerFlow,
    sys::PSY.System;
    timesteps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if timesteps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:timesteps]
        elseif length(timestep_names) != timesteps
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
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    for (bus_no, ix) in bus_lookup
        bus_name = temp_bus_map[bus_no]
        bus = PSY.get_component(PSY.Bus, sys, bus_name)
        bus_type[ix] = PSY.get_bustype(bus)
        if bus_type[ix] == PSY.ACBusTypes.REF
            bus_angles[ix] = 0.0
        else
            bus_angles[ix] = PSY.get_angle(bus)
        end
        bus_magnitude[ix] = PSY.get_magnitude(bus)
    end

    # define injection vectors related to the first timestep
    bus_activepower_injection = zeros(Float64, n_buses)
    bus_reactivepower_injection = zeros(Float64, n_buses)
    _get_injections!(
        bus_activepower_injection,
        bus_reactivepower_injection,
        bus_lookup,
        sys,
    )

    bus_activepower_withdrawals = zeros(Float64, n_buses)
    bus_reactivepower_withdrawals = zeros(Float64, n_buses)
    _get_withdrawals!(
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
    bus_magnitude_1 = zeros(n_buses, 1)
    bus_angles_1 = deepcopy(init_1)
    branch_flow_values_1 = deepcopy(init_2)

    # initial values related to first timestep allocated in the first column
    bus_activepower_injection_1[:, 1] .= bus_activepower_injection
    bus_reactivepower_injection_1[:, 1] .= bus_reactivepower_injection
    bus_activepower_withdrawals_1[:, 1] .= bus_activepower_withdrawals
    bus_reactivepower_withdrawals_1[:, 1] .= bus_reactivepower_withdrawals
    bus_magnitude_1[:, 1] .= bus_magnitude  # for DC case same value accross all timesteps
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        Vector{Float64}(),
        bus_type,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:timesteps], timestep_names)),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
        Vector{Set{Int}}(),
    )
end
