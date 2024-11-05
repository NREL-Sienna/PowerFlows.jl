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
        column number of the previosly mentioned matrices) and their names.
- `valid_ix::Vector{Int}`:
        vector containing the indeces of not slack buses
- `power_network_matrix::M`:
        matrix used for the evaluation of either the power flows or bus angles,
        depending on the method considered.
- `aux_network_matrix::N`:
        matrix used for the evaluation of either the power flows or bus angles,
        depending on the method considered.
- `neighbors::Vector{Set{Int}}`: Vector with the sets of adjacent buses.
"""
struct PowerFlowData{
    M <: PNM.PowerNetworkMatrix,
    N <: Union{PNM.PowerNetworkMatrix, Nothing},
}
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{String, Int}
    bus_activepower_injection::Matrix{Float64}
    bus_reactivepower_injection::Matrix{Float64}
    bus_activepower_withdrawals::Matrix{Float64}
    bus_reactivepower_withdrawals::Matrix{Float64}
    bus_reactivepower_bounds::Vector{Vector{Float64}}
    bus_type::Vector{PSY.ACBusTypes}
    branch_type::Vector{DataType}
    bus_magnitude::Matrix{Float64}
    bus_angles::Matrix{Float64}
    branch_flow_values::Matrix{Float64}
    timestep_map::Dict{Int, String}
    valid_ix::Vector{Int}
    power_network_matrix::M
    aux_network_matrix::N
    neighbors::Vector{Set{Int}}
end

get_bus_lookup(pfd::PowerFlowData) = pfd.bus_lookup
get_branch_lookup(pfd::PowerFlowData) = pfd.branch_lookup
get_bus_activepower_injection(pfd::PowerFlowData) = pfd.bus_activepower_injection
get_bus_reactivepower_injection(pfd::PowerFlowData) = pfd.bus_reactivepower_injection
get_bus_activepower_withdrawals(pfd::PowerFlowData) = pfd.bus_activepower_withdrawals
get_bus_reactivepower_withdrawals(pfd::PowerFlowData) = pfd.bus_reactivepower_withdrawals
get_bus_reactivepower_bounds(pfd::PowerFlowData) = pfd.bus_reactivepower_bounds
get_bus_type(pfd::PowerFlowData) = pfd.bus_type
get_branch_type(pfd::PowerFlowData) = pfd.branch_type
get_bus_magnitude(pfd::PowerFlowData) = pfd.bus_magnitude
get_bus_angles(pfd::PowerFlowData) = pfd.bus_angles
get_branch_flow_values(pfd::PowerFlowData) = pfd.branch_flow_values
get_timestep_map(pfd::PowerFlowData) = pfd.timestep_map
get_valid_ix(pfd::PowerFlowData) = pfd.valid_ix
get_power_network_matrix(pfd::PowerFlowData) = pfd.power_network_matrix
get_aux_network_matrix(pfd::PowerFlowData) = pfd.aux_network_matrix
get_neighbor(pfd::PowerFlowData) = pfd.neighbors

function clear_injection_data!(pfd::PowerFlowData)
    pfd.bus_activepower_injection[:] = 0.0
    pfd.bus_reactivepower_injection[:] = 0.0
    pfd.bus_activepower_withdrawals[:] = 0.0
    pfd.bus_reactivepower_withdrawals[:] = 0.0
    return
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
        use ACPowerFlow() to evaluate the AC PF.
- `sys::PSY.System`:
        container storing the system data to consider in the PowerFlowData
        structure.
- `time_steps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "time_steps". Default
        value = String[].
- `check_connectivity::Bool`:
        Perform connectivity check on the network matrix. Default value = true.

WARNING: functions for the evaluation of the multi-period AC PF still to be implemented.
"""
function PowerFlowData(
    ::ACPowerFlow,
    sys::PSY.System;
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if time_steps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:time_steps]
        elseif length(timestep_names) != time_steps
            error("timestep_names field must have same length as time_steps")
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
    branch_lookup = Dict{String, Int}()
    branch_types = Vector{DataType}(undef, n_branches)
    for (ix, b) in enumerate(branches)
        branch_lookup[PSY.get_name(b)] = ix
        branch_types[ix] = typeof(b)
    end

    # TODO: bus_type might need to also be a Matrix since the type can change for a particular scenario
    bus_type = Vector{PSY.ACBusTypes}(undef, n_buses)
    bus_angles = zeros(Float64, n_buses)
    bus_magnitude = zeros(Float64, n_buses)
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
    )

    _initialize_bus_data!(
        bus_type,
        bus_angles,
        bus_magnitude,
        temp_bus_map,
        bus_lookup,
        sys,
    )

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

    # define fields as matrices whose number of columns is eqault to the number of time_steps
    bus_activepower_injection_1 = zeros(n_buses, time_steps)
    bus_reactivepower_injection_1 = zeros(n_buses, time_steps)
    bus_activepower_withdrawals_1 = zeros(n_buses, time_steps)
    bus_reactivepower_withdrawals_1 = zeros(n_buses, time_steps)
    bus_magnitude_1 = zeros(n_buses, time_steps)
    bus_angles_1 = zeros(n_buses, time_steps)
    branch_flow_values_1 = zeros(n_branches, time_steps)

    bus_reactivepower_bounds = Vector{Vector{Float64}}(undef, n_buses)
    for i in 1:n_buses
        bus_reactivepower_bounds[i] = [0.0, 0.0]
    end
    _get_reactive_power_bound!(bus_reactivepower_bounds, bus_lookup, sys)

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
        bus_reactivepower_bounds,
        bus_type,
        branch_types,
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
- `time_steps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "time_steps". Default
        value = String[].
- `check_connectivity::Bool`:
        Perform connectivity check on the network matrix. Default value = true.
"""
function PowerFlowData(
    ::DCPowerFlow,
    sys::PSY.System;
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if length(timestep_names) == 0
        timestep_names = [string(i) for i in 1:time_steps]
    elseif length(timestep_names) != time_steps
        error("timestep_names field must have same length as time_steps")
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

    branch_types = Vector{DataType}(undef, length(branch_lookup))
    for (ix, b) in enumerate(PNM.get_ac_branches(sys))
        branch_types[ix] = typeof(b)
    end

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
    
    _initialize_bus_data!(
        bus_type,
        bus_angles,
        bus_magnitude,
        temp_bus_map,
        bus_lookup,
        sys,
    )

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
    init_1 = zeros(n_buses, time_steps)
    init_2 = zeros(n_branches, time_steps)

    # define fields as matrices whose number of columns is eqault to the number of time_steps
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
    bus_magnitude_1[:, 1] .= bus_magnitude  # for DC case same value accross all time_steps
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        Vector{Vector{Float64}}(),
        bus_type,
        branch_types,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:time_steps], timestep_names)),
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
- `time_steps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "time_steps". Default
        value = String[].
"""
function PowerFlowData(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if time_steps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:time_steps]
        elseif length(timestep_names) != time_steps
            error("timestep_names field must have same length as time_steps")
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

    branch_types = Vector{DataType}(undef, length(branch_lookup))
    for (ix, b) in enumerate(PNM.get_ac_branches(sys))
        branch_types[ix] = typeof(b)
    end

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

    _initialize_bus_data!(
        bus_type,
        bus_angles,
        bus_magnitude,
        temp_bus_map,
        bus_lookup,
        sys,
    )

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
    init_1 = zeros(n_buses, time_steps)
    init_2 = zeros(n_branches, time_steps)

    # define fields as matrices whose number of columns is eqault to the number of time_steps
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
    bus_magnitude_1[:, 1] .= bus_magnitude  # for DC case same value accross all time_steps
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        Vector{Vector{Float64}}(),
        bus_type,
        branch_types,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:time_steps], timestep_names)),
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
- `time_steps::Int`:
        number of time periods to consider in the PowerFlowData structure. It
        defines the number of columns of the matrices used to store data.
        Default value = 1.
- `timestep_names::Vector{String}`:
        names of the time periods defines by the argmunet "time_steps". Default
        value = String[].
"""
function PowerFlowData(
    ::vPTDFDCPowerFlow,
    sys::PSY.System;
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true)

    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if time_steps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:time_steps]
        elseif length(timestep_names) != time_steps
            error("timestep_names field must have same length as time_steps")
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

    branch_types = Vector{DataType}(undef, length(branch_lookup))
    for (ix, b) in enumerate(PNM.get_ac_branches(sys))
        branch_types[ix] = typeof(b)
    end

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

    _initialize_bus_data!(
        bus_type,
        bus_angles,
        bus_magnitude,
        temp_bus_map,
        bus_lookup,
        sys,
    )

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
    init_1 = zeros(n_buses, time_steps)
    init_2 = zeros(n_branches, time_steps)

    # define fields as matrices whose number of columns is eqault to the number of time_steps
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
    bus_magnitude_1[:, 1] .= bus_magnitude  # for DC case same value accross all time_steps
    bus_angles_1[:, 1] .= bus_angles
    branch_flow_values_1[:, 1] .= zeros(n_branches)

    return PowerFlowData(
        bus_lookup,
        branch_lookup,
        bus_activepower_injection_1,
        bus_reactivepower_injection_1,
        bus_activepower_withdrawals_1,
        bus_reactivepower_withdrawals_1,
        Vector{Vector{Float64}}(),
        bus_type,
        branch_types,
        bus_magnitude_1,
        bus_angles_1,
        branch_flow_values_1,
        Dict(zip([i for i in 1:time_steps], timestep_names)),
        setdiff(1:n_buses, aux_network_matrix.ref_bus_positions),
        power_network_matrix,
        aux_network_matrix,
        Vector{Set{Int}}(),
    )
end
