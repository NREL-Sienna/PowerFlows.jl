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
    power_network_matrix::M
    aux_network_matrix::N
    neighbors::Vector{Set{Int}}
    converged::BitVector
    loss_factors::Union{Matrix{Float64}, Nothing}
    calculate_loss_factors::Bool
    voltage_stability_factors::Union{Matrix{Float64}, Nothing}
    calculate_voltage_stability_factors::Bool
    linear_solver::Symbol  # :klu or :cusolver
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
get_metadata_matrix(pfd::ACPowerFlowData) = pfd.power_network_matrix

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
get_metadata_matrix(pfd::Union{PTDFPowerFlowData, vPTDFPowerFlowData}) =
    pfd.power_network_matrix

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
get_metadata_matrix(pfd::ABAPowerFlowData) = pfd.aux_network_matrix

# true getters for fields:
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
get_bus_reactivepower_bounds(pfd::PowerFlowData) = pfd.bus_reactivepower_bounds
get_bus_slack_participation_factors(pfd::PowerFlowData) =
    pfd.bus_slack_participation_factors
get_generator_slack_participation_factors(pfd::PowerFlowData) =
    pfd.generator_slack_participation_factors
get_bus_type(pfd::PowerFlowData) = pfd.bus_type
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
get_linear_solver(pfd::PowerFlowData) = pfd.linear_solver

# auxiliary getters for the fields of PowerNetworkMatrices we're storing:
# most things we patch through to calls on the metadata matrix:
get_bus_lookup(pfd::PowerFlowData) = PNM.get_bus_lookup(get_metadata_matrix(pfd))
get_bus_axis(pfd::PowerFlowData) = PNM.get_bus_axis(get_metadata_matrix(pfd))
get_arc_lookup(pfd::PowerFlowData) = PNM.get_arc_lookup(get_metadata_matrix(pfd))
get_arc_axis(pfd::PowerFlowData) = PNM.get_arc_axis(get_metadata_matrix(pfd))
get_network_reduction_data(pfd::PowerFlowData) =
    PNM.get_network_reduction_data(get_metadata_matrix(pfd))
get_valid_ix(pdf::PowerFlowData) = Not(PNM.get_ref_bus_position(get_metadata_matrix(pdf)))

# the ybus matrix itself doesn't have an arc axis, so we have to special-case it.
get_arc_lookup(pfd::ACPowerFlowData) =
    PNM.get_arc_lookup(pfd.power_network_matrix.arc_admittance_from_to)
get_arc_axis(pfd::ACPowerFlowData) =
    PNM.get_arc_axis(pfd.power_network_matrix.arc_admittance_from_to)

# so we can initialize things to the correct size inside the below constructor.
# No `PowerFlowData` instance, so can't call get_arc_axis or similar to get the size.
arc_count(::ACPowerFlow,
    power_network_matrix::PNM.PowerNetworkMatrix,
    ::Union{PNM.PowerNetworkMatrix, Nothing}) = length(PNM.get_arc_axis(power_network_matrix.arc_admittance_from_to))
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

"""
Sets the two `PowerNetworkMatrix` fields and a few others (`timesteps`, `timestep_map`), 
then creates arrays of default values (usually zeros) for the rest.
"""
function PowerFlowData(
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
    linear_solver = get_linear_solver(pf)
    if !isnothing(get_slack_participation_factors(pf))
        empty_slack_participation_factors = Dict{Tuple{DataType, String}, Float64}[]
        sizehint!(empty_slack_participation_factors, n_timesteps)
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
        linear_solver,
    )
end

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

# NOTE: remove this once network reductions are fully implemented
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

function make_and_initialize_powerflow_data(
    pf::PowerFlowEvaluationModel,
    sys::PSY.System,
    power_network_matrix::M,
    aux_network_matrix::N;
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    neighbors = Vector{Set{Int}}(),
    correct_bustypes::Bool = false,
) where {M <: PNM.PowerNetworkMatrix, N <: Union{PNM.PowerNetworkMatrix, Nothing}}
    data = PowerFlowData(
        pf,
        power_network_matrix,
        aux_network_matrix,
        time_steps;
        timestep_names = timestep_names,
        neighbors = neighbors,
    )
    initialize_powerflow_data!(data, pf, sys; correct_bustypes = correct_bustypes)
    return data
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

WARNING: functions for the evaluation of the multi-period AC PF still to be implemented.
"""
function PowerFlowData(
    pf::ACPowerFlow{<:ACPowerFlowSolverType},
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes::Bool = false,
)
    network_reduction_message(network_reductions, pf)
    power_network_matrix = PNM.Ybus(
        sys;
        network_reductions = network_reductions,
        make_arc_admittance_matrices = true,
        include_constant_impedance_loads = false,
    )
    neighbors = _calculate_neighbors(power_network_matrix)

    if get_robust_power_flow(pf)
        aux_network_matrix = PNM.ABA_Matrix(sys; network_reductions = network_reductions)
    else
        aux_network_matrix = nothing
    end

    return make_and_initialize_powerflow_data(
        pf,
        sys,
        power_network_matrix,
        aux_network_matrix;
        time_steps = time_steps,
        timestep_names = timestep_names,
        neighbors = neighbors,
        correct_bustypes = correct_bustypes,
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
"""
function PowerFlowData(
    pf::DCPowerFlow,
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes = false,
)
    network_reduction_message(network_reductions, DCPowerFlow())
    # get the network matrices
    power_network_matrix =
        PNM.ABA_Matrix(sys; factorize = true, network_reductions = network_reductions)
    aux_network_matrix = PNM.BA_Matrix(sys; network_reductions = network_reductions)
    return make_and_initialize_powerflow_data(
        pf,
        sys,
        power_network_matrix,
        aux_network_matrix;
        time_steps = time_steps,
        timestep_names = timestep_names,
        correct_bustypes = correct_bustypes,
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
    pf::PTDFDCPowerFlow,
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes = false,
)
    network_reduction_message(network_reductions, PTDFDCPowerFlow())
    # get the network matrices
    power_network_matrix = PNM.PTDF(sys; network_reductions = network_reductions)
    aux_network_matrix =
        PNM.ABA_Matrix(sys; factorize = true, network_reductions = network_reductions)
    return make_and_initialize_powerflow_data(
        pf,
        sys,
        power_network_matrix,
        aux_network_matrix;
        time_steps = time_steps,
        timestep_names = timestep_names,
        correct_bustypes = correct_bustypes,
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
    pf::vPTDFDCPowerFlow,
    sys::PSY.System;
    network_reductions::Vector{PNM.NetworkReduction} = Vector{PNM.NetworkReduction}(),
    time_steps::Int = 1,
    timestep_names::Vector{String} = String[],
    correct_bustypes = false)
    network_reduction_message(network_reductions, vPTDFDCPowerFlow())

    # get the network matrices
    power_network_matrix = PNM.VirtualPTDF(sys; network_reductions = network_reductions) # evaluates an empty virtual PTDF
    aux_network_matrix =
        PNM.ABA_Matrix(sys; factorize = true, network_reductions = network_reductions)

    return make_and_initialize_powerflow_data(
        pf,
        sys,
        power_network_matrix,
        aux_network_matrix;
        time_steps = time_steps,
        timestep_names = timestep_names,
        correct_bustypes = correct_bustypes,
    )
end

"""
Create an appropriate `PowerFlowContainer` for the given `PowerFlowEvaluationModel` and initialize it from the given `PSY.System`.

# Arguments:
- `pfem::PowerFlowEvaluationModel`: power flow model to construct a container for (e.g., `DCPowerFlow()`)
- `sys::PSY.System`: the system from which to initialize the power flow container
- `time_steps::Int`: number of time periods to consider (default is `1`)
- `timestep_names::Vector{String}`: names of the time periods defines by the argument "time_steps". Default value is `String[]`.
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
