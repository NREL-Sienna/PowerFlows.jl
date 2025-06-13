abstract type PowerFlowContainer end

"""
Trait signifying whether the `PowerFlowContainer` can represent multi-period data. Must be
implemented for all concrete subtypes.
"""
supports_multi_period(x::PowerFlowContainer) =
    throw(
        IS.NotImplementedError(
            "supports_multi_period must be implemented for $(typeof(x))"),
    )

"A `PowerFlowContainer` that represents its data as a `PSY.System`"
abstract type SystemPowerFlowContainer <: PowerFlowContainer end

get_system(container::SystemPowerFlowContainer) = container.system

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
- `bus_reactivepower_bounds::Matrix{Float64}`:
        "(b, t)" matrix containing upper and lower bounds for the reactive supply at each
        bus at each time period.
- `bus_type::Matrix{PSY.ACBusTypes}`:
        "(b, t)" matrix containing type of buses present in the system, ordered
        according to "bus_lookup," at each time period.
- `bus_magnitude::Matrix{Float64}`:
        "(b, t)" matrix containing the bus magnitudes, ordered according to
        "bus_lookup". b: number of buses, t: number of time period.
- `bus_angles::Matrix{Float64}`:
        "(b, t)" matrix containing the bus angles, ordered according to
        "bus_lookup". b: number of buses, t: number of time period.
- `branch_activepower_flow_from_to::Matrix{Float64}`:
        "(br, t)" matrix containing the active power flows measured at the `from` bus,
        ordered according to "branch_lookup". br: number of branches, t: number of time
        period.
- `branch_reactivepower_flow_from_to::Matrix{Float64}`:
        "(br, t)" matrix containing the reactive power flows measured at the `from` bus,
        ordered according to "branch_lookup". br: number of branches, t: number of time
        period.
- `branch_activepower_flow_to_from::Matrix{Float64}`:
        "(br, t)" matrix containing the active power flows measured at the `to` bus, ordered
        according to "branch_lookup". br: number of branches, t: number of time period.
- `branch_reactivepower_flow_to_from::Matrix{Float64}`:
        "(br, t)" matrix containing the reactive power flows measured at the `to` bus,
        ordered according to "branch_lookup". br: number of branches, t: number of time
        period.
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
} <: PowerFlowContainer
    bus_lookup::Dict{Int, Int}
    branch_lookup::Dict{String, Int}
    bus_activepower_injection::Matrix{Float64}
    bus_reactivepower_injection::Matrix{Float64}
    bus_activepower_withdrawals::Matrix{Float64}
    bus_reactivepower_withdrawals::Matrix{Float64}
    bus_reactivepower_bounds::Matrix{Vector{Float64}}
    generator_slack_participation_factors::Union{
        Vector{Dict{Tuple{DataType, String}, Float64}},
        Nothing,
    }
    bus_slack_participation_factors::SparseMatrixCSC{Float64, Int}
    bus_type::Matrix{PSY.ACBusTypes}
    branch_type::Vector{DataType}
    bus_magnitude::Matrix{Float64}
    bus_angles::Matrix{Float64}
    branch_activepower_flow_from_to::Matrix{Float64}
    branch_reactivepower_flow_from_to::Matrix{Float64}
    branch_activepower_flow_to_from::Matrix{Float64}
    branch_reactivepower_flow_to_from::Matrix{Float64}
    timestep_map::Dict{Int, String}
    valid_ix::Vector{Int}
    power_network_matrix::M
    aux_network_matrix::N
    neighbors::Vector{Set{Int}}
    converged::Vector{Bool}
    loss_factors::Union{Matrix{Float64}, Nothing}
    calculate_loss_factors::Bool
end

get_bus_lookup(pfd::PowerFlowData) = pfd.bus_lookup
get_branch_lookup(pfd::PowerFlowData) = pfd.branch_lookup
get_bus_activepower_injection(pfd::PowerFlowData) = pfd.bus_activepower_injection
get_bus_reactivepower_injection(pfd::PowerFlowData) = pfd.bus_reactivepower_injection
get_bus_activepower_withdrawals(pfd::PowerFlowData) = pfd.bus_activepower_withdrawals
get_bus_reactivepower_withdrawals(pfd::PowerFlowData) = pfd.bus_reactivepower_withdrawals
get_bus_reactivepower_bounds(pfd::PowerFlowData) = pfd.bus_reactivepower_bounds
get_bus_slack_participation_factors(pfd::PowerFlowData) =
    pfd.bus_slack_participation_factors
get_generator_slack_participation_factors(pfd::PowerFlowData) =
    pfd.generator_slack_participation_factors
get_bus_type(pfd::PowerFlowData) = pfd.bus_type
get_branch_type(pfd::PowerFlowData) = pfd.branch_type
get_bus_magnitude(pfd::PowerFlowData) = pfd.bus_magnitude
get_bus_angles(pfd::PowerFlowData) = pfd.bus_angles
get_branch_activepower_flow_from_to(pfd::PowerFlowData) =
    pfd.branch_activepower_flow_from_to
get_branch_reactivepower_flow_from_to(pfd::PowerFlowData) =
    pfd.branch_reactivepower_flow_from_to
get_branch_activepower_flow_to_from(pfd::PowerFlowData) =
    pfd.branch_activepower_flow_to_from
get_branch_reactivepower_flow_to_from(pfd::PowerFlowData) =
    pfd.branch_reactivepower_flow_to_from
get_timestep_map(pfd::PowerFlowData) = pfd.timestep_map
get_valid_ix(pfd::PowerFlowData) = pfd.valid_ix
get_power_network_matrix(pfd::PowerFlowData) = pfd.power_network_matrix
get_aux_network_matrix(pfd::PowerFlowData) = pfd.aux_network_matrix
get_neighbor(pfd::PowerFlowData) = pfd.neighbors
supports_multi_period(::PowerFlowData) = true
get_converged(pfd::PowerFlowData) = pfd.converged
get_loss_factors(pfd::PowerFlowData) = pfd.loss_factors

function clear_injection_data!(pfd::PowerFlowData)
    pfd.bus_activepower_injection .= 0.0
    pfd.bus_reactivepower_injection .= 0.0
    pfd.bus_activepower_withdrawals .= 0.0
    pfd.bus_reactivepower_withdrawals .= 0.0
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
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    sys::PSY.System;
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    check_connectivity::Bool = true,
    correct_bustypes = false)
    calculate_loss_factors = pf.calculate_loss_factors
    generator_slack_participation_factors = pf.generator_slack_participation_factors
    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if time_steps != 0
        if length(timestep_names) == 0
            timestep_names = [string(i) for i in 1:time_steps]
        elseif length(timestep_names) != time_steps
            error("timestep_names field must have same length as time_steps")
        end
    end

    timestep_map = Dict(zip([i for i in 1:time_steps], timestep_names))

    # get data for calculations
    power_network_matrix = PNM.Ybus(
        sys;
        check_connectivity = check_connectivity,
        make_branch_admittance_matrices = true,
    )

    # get number of buses and branches
    n_buses = length(axes(power_network_matrix, 1))
    buses = PNM.get_buses(sys)
    ref_bus_positions = PNM.find_slack_positions(buses)

    branches = PNM.get_ac_branches(sys)
    n_branches = length(branches)

    bus_lookup = power_network_matrix.lookup[2]
    branch_lookup = Dict{String, Int}()
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )
    branch_types = Vector{DataType}(undef, n_branches)
    for (ix, b) in enumerate(branches)
        branch_lookup[PSY.get_name(b)] = ix
        branch_types[ix] = typeof(b)
    end

    valid_ix = setdiff(1:n_buses, ref_bus_positions)
    neighbors = _calculate_neighbors(power_network_matrix)
    aux_network_matrix = nothing
    converged = fill(false, time_steps)
    # loss factors order matches the order of buses in the grid model, and is calculated for all buses including ref buses (equals 0 for ref buses)
    if calculate_loss_factors
        loss_factors = Matrix{Float64}(undef, (n_buses, length(timestep_names)))
    else
        loss_factors = nothing
    end

    return make_powerflowdata(
        sys,
        time_steps,
        power_network_matrix,
        aux_network_matrix,
        n_buses,
        n_branches,
        bus_lookup,
        branch_lookup,
        temp_bus_map,
        branch_types,
        timestep_map,
        valid_ix,
        neighbors,
        converged,
        loss_factors,
        calculate_loss_factors,
        correct_bustypes,
        generator_slack_participation_factors,
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
    check_connectivity::Bool = true,
    correct_bustypes = false)

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
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )
    valid_ix = setdiff(1:n_buses, aux_network_matrix.ref_bus_positions)
    converged = fill(false, time_steps)
    loss_factors = nothing
    calculate_loss_factors = false
    return make_dc_powerflowdata(
        sys,
        time_steps,
        timestep_names,
        power_network_matrix,
        aux_network_matrix,
        n_buses,
        n_branches,
        bus_lookup,
        branch_lookup,
        temp_bus_map,
        valid_ix,
        converged,
        loss_factors,
        correct_bustypes,
        calculate_loss_factors,
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
    check_connectivity::Bool = true,
    correct_bustypes = false)

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
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )
    valid_ix = setdiff(1:n_buses, aux_network_matrix.ref_bus_positions)
    converged = fill(false, time_steps)
    loss_factors = nothing
    calculate_loss_factors = false
    return make_dc_powerflowdata(
        sys,
        time_steps,
        timestep_names,
        power_network_matrix,
        aux_network_matrix,
        n_buses,
        n_branches,
        bus_lookup,
        branch_lookup,
        temp_bus_map,
        valid_ix,
        converged,
        loss_factors,
        correct_bustypes,
        calculate_loss_factors,
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
    check_connectivity::Bool = true,
    correct_bustypes = false)

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
    temp_bus_map = Dict{Int, String}(
        PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.ACBus, sys)
    )
    valid_ix = setdiff(1:n_buses, aux_network_matrix.ref_bus_positions)
    converged = fill(false, time_steps)
    loss_factors = nothing
    calculate_loss_factors = false
    return make_dc_powerflowdata(
        sys,
        time_steps,
        timestep_names,
        power_network_matrix,
        aux_network_matrix,
        n_buses,
        n_branches,
        bus_lookup,
        branch_lookup,
        temp_bus_map,
        valid_ix,
        converged,
        loss_factors,
        correct_bustypes,
        calculate_loss_factors,
    )
end

"""
Create an appropriate `PowerFlowContainer` for the given `PowerFlowEvaluationModel` and initialize it from the given `PSY.System`.

# Arguments:
- `pfem::PowerFlowEvaluationModel`: power flow model to construct a container for (e.g., `DCPowerFlow()`)
- `sys::PSY.System`: the system from which to initialize the power flow container
- `time_steps::Int`: number of time periods to consider (default is `1`)
- `timestep_names::Vector{String}`: names of the time periods defines by the argument "time_steps". Default value is `String[]`.
- `check_connectivity::Bool`: Perform connectivity check on the network matrix. Default value is `true`.
"""
function make_power_flow_container end

make_power_flow_container(
    pfem::ACPowerFlow{<:ACPowerFlowSolverType},
    sys::PSY.System;
    kwargs...,
) =
    PowerFlowData(pfem, sys; kwargs...)

make_power_flow_container(pfem::DCPowerFlow, sys::PSY.System; kwargs...) =
    PowerFlowData(pfem, sys; kwargs...)

make_power_flow_container(pfem::PTDFDCPowerFlow, sys::PSY.System; kwargs...) =
    PowerFlowData(pfem, sys; kwargs...)

make_power_flow_container(pfem::vPTDFDCPowerFlow, sys::PSY.System; kwargs...) =
    PowerFlowData(pfem, sys; kwargs...)
