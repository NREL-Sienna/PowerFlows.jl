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

In the below descriptions, "number of buses" should be understood as "number of buses remaining,
after the network reduction." Similarly, we use "arcs" instead of "branches" to distinguish 
between network elements (post-reduction) and system objects (pre-reduction).
# Arguments:
- `bus_lookup::Dict{Int, Int}`:
        dictionary linking the system's bus number with the rows of either
        "power_network_matrix" or "aux_network_matrix".
- `arc_lookup::Dict{Tuple{Int, Int}, Int}`:
        dictionary linking the arc name with the column name of either the
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
- `bus_activepower_constant_current_withdrawals::Matrix{Float64}`:
        "(b, t)" matrix containing the bus active power constant current
        withdrawals. b: number of buses, t: number of time period.
- `bus_reactivepower_constant_current_withdrawals::Matrix{Float64}`:
        "(b, t)" matrix containing the bus reactive power constant current
        withdrawals. b: number of buses, t: number of time period.
- `bus_activepower_constant_impedance_withdrawals::Matrix{Float64}`:
        "(b, t)" matrix containing the bus active power constant impedance
        withdrawals. b: number of buses, t: number of time period.
- `bus_reactivepower_constant_impedance_withdrawals::Matrix{Float64}`:  
        "(b, t)" matrix containing the bus reactive power constant impedance
        withdrawals. b: number of buses, t: number of time period.
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
- `arc_activepower_flow_from_to::Matrix{Float64}`:
        "(br, t)" matrix containing the active power flows measured at the `from` bus,
        ordered according to "arc_lookup". br: number of arcs, t: number of time
        period.
- `arc_reactivepower_flow_from_to::Matrix{Float64}`:
        "(br, t)" matrix containing the reactive power flows measured at the `from` bus,
        ordered according to "arc_lookup". br: number of arcs, t: number of time
        period.
- `arc_activepower_flow_to_from::Matrix{Float64}`:
        "(br, t)" matrix containing the active power flows measured at the `to` bus, ordered
        according to "arc_lookup". br: number of arcs, t: number of time period.
- `arc_reactivepower_flow_to_from::Matrix{Float64}`:
        "(br, t)" matrix containing the reactive power flows measured at the `to` bus,
        ordered according to "arc_lookup". br: number of arcs, t: number of time
        period.
- `timestep_map::Dict{Int, S}`:
        dictionary mapping the number of the time periods (corresponding to the
        column number of the previously mentioned matrices) and their names.
- `valid_ix::Vector{Int}`:
        vector containing the indices of not slack buses
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
    arc_lookup::Dict{Tuple{Int, Int}, Int}
    bus_activepower_injection::Matrix{Float64}
    bus_reactivepower_injection::Matrix{Float64}
    bus_activepower_withdrawals::Matrix{Float64}
    bus_reactivepower_withdrawals::Matrix{Float64}
    bus_activepower_constant_current_withdrawals::Matrix{Float64}
    bus_reactivepower_constant_current_withdrawals::Matrix{Float64}
    bus_activepower_constant_impedance_withdrawals::Matrix{Float64}
    bus_reactivepower_constant_impedance_withdrawals::Matrix{Float64}
    bus_reactivepower_bounds::Matrix{Tuple{Float64, Float64}}
    generator_slack_participation_factors::Union{
        Vector{Dict{Tuple{DataType, String}, Float64}},
        Nothing,
    }
    bus_slack_participation_factors::SparseMatrixCSC{Float64, Int}
    bus_type::Matrix{PSY.ACBusTypes}
    bus_magnitude::Matrix{Float64}
    bus_angles::Matrix{Float64}
    arc_activepower_flow_from_to::Matrix{Float64}
    arc_reactivepower_flow_from_to::Matrix{Float64}
    arc_activepower_flow_to_from::Matrix{Float64}
    arc_reactivepower_flow_to_from::Matrix{Float64}
    timestep_map::Dict{Int, String}
    valid_ix::Vector{Int}
    power_network_matrix::M
    aux_network_matrix::N
    neighbors::Vector{Set{Int}}
    converged::Vector{Bool}
    loss_factors::Union{Matrix{Float64}, Nothing}
    calculate_loss_factors::Bool
    voltage_stability_factors::Union{Matrix{Float64}, Nothing}
    calculate_voltage_stability_factors::Bool
end

# aliases for specific type parameter combinations.
const ACPowerFlowData = PowerFlowData{
    PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
    <:Union{
        PNM.ABA_Matrix{Tuple{Vector{Int64}, Vector{Int64}},
            Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
            Nothing},
        Nothing,
    },
}
get_arc_axis(pfd::ACPowerFlowData) =
    PNM.get_arc_axis(pfd.power_network_matrix.branch_admittance_from_to)
get_bus_axis(pfd::ACPowerFlowData) =
    PNM.get_bus_axis(pfd.power_network_matrix)
get_network_reduction_data(pfd::ACPowerFlowData) =
    PNM.get_network_reduction_data(pfd.power_network_matrix)

const PTDFPowerFlowData = PowerFlowData{
    PNM.PTDF{
        Tuple{Vector{Int64}, Vector{Tuple{Int, Int}}},
        Tuple{Dict{Int64, Int64}, Dict{Tuple{Int, Int}, Int64}},
        Matrix{Float64},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

const vPTDFPowerFlowData = PowerFlowData{
    PNM.VirtualPTDF{
        Tuple{Vector{Tuple{Int, Int}}, Vector{Int64}},
        Tuple{Dict{Tuple{Int, Int}, Int64}, Dict{Int64, Int64}},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}
get_arc_axis(pfd::Union{PTDFPowerFlowData, vPTDFPowerFlowData}) =
    PNM.get_arc_axis(pfd.power_network_matrix)
get_bus_axis(pfd::Union{PTDFPowerFlowData, vPTDFPowerFlowData}) =
    PNM.get_bus_axis(pfd.power_network_matrix)
get_network_reduction_data(pfd::Union{PTDFPowerFlowData, vPTDFPowerFlowData}) =
    PNM.get_network_reduction_data(pfd.power_network_matrix)

const ABAPowerFlowData = PowerFlowData{
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
    PNM.BA_Matrix{
        Tuple{Vector{Int64}, Vector{Tuple{Int, Int}}},
        Tuple{Dict{Int64, Int64}, Dict{Tuple{Int, Int}, Int64}}},
}
get_arc_axis(pfd::ABAPowerFlowData) = PNM.get_arc_axis(pfd.aux_network_matrix)
get_bus_axis(pfd::ABAPowerFlowData) = PNM.get_bus_axis(pfd.aux_network_matrix)
get_network_reduction_data(pfd::ABAPowerFlowData) =
    PNM.get_network_reduction_data(pfd.aux_network_matrix)

get_bus_lookup(pfd::PowerFlowData) = pfd.bus_lookup
get_arc_lookup(pfd::PowerFlowData) = pfd.arc_lookup
get_bus_activepower_injection(pfd::PowerFlowData) = pfd.bus_activepower_injection
get_bus_reactivepower_injection(pfd::PowerFlowData) = pfd.bus_reactivepower_injection
get_bus_activepower_withdrawals(pfd::PowerFlowData) = pfd.bus_activepower_withdrawals
get_bus_activepower_constant_current_withdrawals(pfd::PowerFlowData) =
    pfd.bus_activepower_constant_current_withdrawals
get_bus_activepower_constant_impedance_withdrawals(pfd::PowerFlowData) =
    pfd.bus_activepower_constant_impedance_withdrawals
get_bus_reactivepower_withdrawals(pfd::PowerFlowData) = pfd.bus_reactivepower_withdrawals
get_bus_reactivepower_constant_current_withdrawals(pfd::PowerFlowData) =
    pfd.bus_reactivepower_constant_current_withdrawals
get_bus_reactivepower_constant_impedance_withdrawals(pfd::PowerFlowData) =
    pfd.bus_reactivepower_constant_impedance_withdrawals

function get_bus_activepower_total_withdrawals(pfd::PowerFlowData, ix::Int, time_step::Int)
    return pfd.bus_activepower_withdrawals[ix, time_step] +
           pfd.bus_activepower_constant_current_withdrawals[ix, time_step] *
           pfd.bus_magnitude[ix, time_step] +
           pfd.bus_activepower_constant_impedance_withdrawals[ix, time_step] *
           pfd.bus_magnitude[ix, time_step]^2
end

function get_bus_reactivepower_total_withdrawals(
    pfd::PowerFlowData,
    ix::Int,
    time_step::Int,
)
    return pfd.bus_reactivepower_withdrawals[ix, time_step] +
           pfd.bus_reactivepower_constant_current_withdrawals[ix, time_step] *
           pfd.bus_magnitude[ix, time_step] +
           pfd.bus_reactivepower_constant_impedance_withdrawals[ix, time_step] *
           pfd.bus_magnitude[ix, time_step]^2
end

get_bus_reactivepower_bounds(pfd::PowerFlowData) = pfd.bus_reactivepower_bounds
get_bus_slack_participation_factors(pfd::PowerFlowData) =
    pfd.bus_slack_participation_factors
get_generator_slack_participation_factors(pfd::PowerFlowData) =
    pfd.generator_slack_participation_factors
get_bus_type(pfd::PowerFlowData) = pfd.bus_type
# get_arc_type(pfd::PowerFlowData) = pfd.branch_type
get_bus_magnitude(pfd::PowerFlowData) = pfd.bus_magnitude
get_bus_angles(pfd::PowerFlowData) = pfd.bus_angles
get_arc_activepower_flow_from_to(pfd::PowerFlowData) =
    pfd.arc_activepower_flow_from_to
get_arc_reactivepower_flow_from_to(pfd::PowerFlowData) =
    pfd.arc_reactivepower_flow_from_to
get_arc_activepower_flow_to_from(pfd::PowerFlowData) =
    pfd.arc_activepower_flow_to_from
get_arc_reactivepower_flow_to_from(pfd::PowerFlowData) =
    pfd.arc_reactivepower_flow_to_from
get_timestep_map(pfd::PowerFlowData) = pfd.timestep_map
get_valid_ix(pfd::PowerFlowData) = pfd.valid_ix
get_power_network_matrix(pfd::PowerFlowData) = pfd.power_network_matrix
get_aux_network_matrix(pfd::PowerFlowData) = pfd.aux_network_matrix
get_neighbor(pfd::PowerFlowData) = pfd.neighbors
supports_multi_period(::PowerFlowData) = true
get_converged(pfd::PowerFlowData) = pfd.converged
get_loss_factors(pfd::PowerFlowData) = pfd.loss_factors
get_calculate_loss_factors(pfd::PowerFlowData) = pfd.calculate_loss_factors
get_voltage_stability_factors(pfd::PowerFlowData) = pfd.voltage_stability_factors
get_calculate_voltage_stability_factors(pfd::PowerFlowData) =
    pfd.calculate_voltage_stability_factors

function clear_injection_data!(pfd::PowerFlowData)
    pfd.bus_activepower_injection .= 0.0
    pfd.bus_reactivepower_injection .= 0.0
    pfd.bus_activepower_withdrawals .= 0.0
    pfd.bus_activepower_constant_current_withdrawals .= 0.0
    pfd.bus_activepower_constant_impedance_withdrawals .= 0.0
    pfd.bus_reactivepower_withdrawals .= 0.0
    pfd.bus_reactivepower_constant_current_withdrawals .= 0.0
    pfd.bus_reactivepower_constant_impedance_withdrawals .= 0.0
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

function network_reduction_message(
    nrs::Vector{PNM.NetworkReduction},
    m::PowerFlowEvaluationModel,
)
    if any(isa.(nrs, (PNM.WardReduction,)))
        throw(IS.NotImplementedError("Ward reduction is not supported yet."))
    end
    if m isa ACPowerFlow && any(isa.(nrs, (PNM.RadialReduction,)))
        @error "AC Power Flow with Radial Network Reduction: feature is a work-in-progress. The power flow will likely fail to converge."
    end
    if any(isa.(nrs, (PNM.DegreeTwoReduction,)))
        @warn "Degree 2 network reductions mis-report branch power flows, but bus voltage results are correct. Use with caution."
    end
    return
end

"""
Function for the definition of the PowerFlowData structure given the System
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
        names of the time periods defines by the argument "time_steps". Default
        value = String[].
- `check_connectivity::Bool`:
        Perform connectivity check on the network matrix. Default value = true.

WARNING: functions for the evaluation of the multi-period AC PF still to be implemented.
"""
function PowerFlowData(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes::Bool = false)
    network_reduction_message(network_reductions, pf)
    calculate_loss_factors = pf.calculate_loss_factors
    generator_slack_participation_factors = pf.generator_slack_participation_factors
    calculate_voltage_stability_factors = pf.calculate_voltage_stability_factors
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
        network_reductions = network_reductions,
        make_branch_admittance_matrices = true,
        include_constant_impedance_loads = false,
    )

    # get number of arcs and branches
    arc_lookup = PNM.get_arc_lookup(power_network_matrix.branch_admittance_from_to)
    bus_lookup = PNM.get_bus_lookup(power_network_matrix)
    n_buses = length(bus_lookup)

    ref_bus_positions = PNM.get_ref_bus_position(power_network_matrix)

    valid_ix = setdiff(1:n_buses, ref_bus_positions)
    neighbors = _calculate_neighbors(power_network_matrix)
    converged = fill(false, time_steps)
    # PERF: type instability.
    # loss factors order matches the order of buses in the grid model, and is calculated for all buses including ref buses (equals 0 for ref buses)
    if calculate_loss_factors
        loss_factors = Matrix{Float64}(undef, (n_buses, length(timestep_names)))
    else
        loss_factors = nothing
    end
    if calculate_voltage_stability_factors
        voltage_stability_factors =
            Matrix{Float64}(undef, (n_buses, length(timestep_names)))
    else
        voltage_stability_factors = nothing
    end
    if get_robust_power_flow(pf)
        aux_network_matrix = PNM.ABA_Matrix(sys)
    else
        aux_network_matrix = nothing
    end
    return make_powerflowdata(
        sys,
        time_steps,
        power_network_matrix,
        aux_network_matrix,
        bus_lookup,
        arc_lookup,
        timestep_map,
        valid_ix,
        neighbors,
        converged,
        loss_factors,
        calculate_loss_factors,
        voltage_stability_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        correct_bustypes,
    )
end

# DC Power Flow Data based on ABA and BA matrices
"""
Function for the definition of the PowerFlowData structure given the System
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
        names of the time periods defines by the argument "time_steps". Default
        value = String[].
- `check_connectivity::Bool`:
        Perform connectivity check on the network matrix. Default value = true.
"""
function PowerFlowData(
    ::DCPowerFlow,
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes = false)
    network_reduction_message(network_reductions, DCPowerFlow())
    # assign timestep_names
    # timestep names are then allocated in a dictionary to map matrix columns
    if length(timestep_names) == 0
        timestep_names = [string(i) for i in 1:time_steps]
    elseif length(timestep_names) != time_steps
        error("timestep_names field must have same length as time_steps")
    end

    # get the network matrices
    power_network_matrix =
        PNM.ABA_Matrix(sys; factorize = true, network_reductions = network_reductions)
    aux_network_matrix = PNM.BA_Matrix(sys; network_reductions = network_reductions)

    # get number of arcs and branches
    bus_lookup = PNM.get_bus_lookup(aux_network_matrix)
    arc_lookup = PNM.get_arc_lookup(aux_network_matrix)

    valid_ix = setdiff(1:length(bus_lookup), PNM.get_ref_bus_position(aux_network_matrix))
    converged = fill(false, time_steps)
    loss_factors = nothing
    calculate_loss_factors = false
    generator_slack_participation_factors = nothing
    voltage_stability_factors = nothing
    calculate_voltage_stability_factors = false
    return make_dc_powerflowdata(
        sys,
        time_steps,
        timestep_names,
        power_network_matrix,
        aux_network_matrix,
        bus_lookup,
        arc_lookup,
        valid_ix,
        converged,
        loss_factors,
        calculate_loss_factors,
        voltage_stability_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        correct_bustypes,
    )
end

# DC Power Flow Data with PTDF matrix
"""
Function for the definition of the PowerFlowData structure given the System
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
        names of the time periods defines by the argument "time_steps". Default
        value = String[].
"""

function PowerFlowData(
    ::PTDFDCPowerFlow,
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes = false)
    network_reduction_message(network_reductions, PTDFDCPowerFlow())
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
    power_network_matrix = PNM.PTDF(sys; network_reductions = network_reductions)
    aux_network_matrix =
        PNM.ABA_Matrix(sys; factorize = true, network_reductions = network_reductions)

    # get number of arcs and branches
    bus_lookup = PNM.get_bus_lookup(power_network_matrix)
    arc_lookup = PNM.get_arc_lookup(power_network_matrix)
    valid_ix = setdiff(1:length(bus_lookup), PNM.get_ref_bus_position(power_network_matrix))
    converged = fill(false, time_steps)
    loss_factors = nothing
    calculate_loss_factors = false
    generator_slack_participation_factors = nothing
    voltage_stability_factors = nothing
    calculate_voltage_stability_factors = false
    return make_dc_powerflowdata(
        sys,
        time_steps,
        timestep_names,
        power_network_matrix,
        aux_network_matrix,
        bus_lookup,
        arc_lookup,
        valid_ix,
        converged,
        loss_factors,
        calculate_loss_factors,
        voltage_stability_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        correct_bustypes,
    )
end

# DC Power Flow Data with virtual PTDF matrix
"""
Function for the definition of the PowerFlowData structure given the System
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
        names of the time periods defines by the argument "time_steps". Default
        value = String[].
"""
function PowerFlowData(
    ::vPTDFDCPowerFlow,
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes = false)
    network_reduction_message(network_reductions, vPTDFDCPowerFlow())
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
    power_network_matrix = PNM.VirtualPTDF(sys; network_reductions = network_reductions) # evaluates an empty virtual PTDF
    aux_network_matrix =
        PNM.ABA_Matrix(sys; factorize = true, network_reductions = network_reductions)

    bus_lookup = power_network_matrix.lookup[2]
    arc_lookup = power_network_matrix.lookup[1]
    valid_ix = setdiff(1:length(bus_lookup), PNM.get_ref_bus_position(power_network_matrix))
    converged = fill(false, time_steps)
    loss_factors = nothing
    calculate_loss_factors = false
    generator_slack_participation_factors = nothing
    voltage_stability_factors = nothing
    calculate_voltage_stability_factors = false
    return make_dc_powerflowdata(
        sys,
        time_steps,
        timestep_names,
        power_network_matrix,
        aux_network_matrix,
        bus_lookup,
        arc_lookup,
        valid_ix,
        converged,
        loss_factors,
        calculate_loss_factors,
        voltage_stability_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        correct_bustypes,
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
