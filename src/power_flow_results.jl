"""
Abstract supertype for all power flow result containers.
"""
abstract type AbstractPowerFlowResults end

"""Fields stored in result containers and forwarded from `PowerFlowData` via `getproperty`."""
const _RESULT_FIELD_NAMES = (
    :bus_magnitude,
    :bus_angles,
    :bus_type,
    :bus_active_power_injections,
    :bus_reactive_power_injections,
    :bus_active_power_withdrawals,
    :bus_reactive_power_withdrawals,
    :arc_active_power_flow_from_to,
    :arc_reactive_power_flow_from_to,
    :arc_active_power_flow_to_from,
    :arc_reactive_power_flow_to_from,
    :arc_angle_differences,
    :converged,
    :loss_factors,
    :voltage_stability_factors,
    :arc_active_power_losses,
    :lcc_active_power_flow_from_to,
    :lcc_active_power_flow_to_from,
    :lcc_reactive_power_flow_from_to,
    :lcc_reactive_power_flow_to_from,
)

"""
    TimePowerFlowData <: AbstractPowerFlowResults

Container for time-series power flow results, storing bus and arc quantities across
multiple time steps as matrices (rows = buses/arcs, columns = time steps).

# Fields
- `bus_magnitude::Matrix{Float64}`: bus voltage magnitudes.
- `bus_angles::Matrix{Float64}`: bus voltage angles.
- `bus_type::Matrix{PSY.ACBusTypes}`: bus type classification at each time step.
- `bus_active_power_injections::Matrix{Float64}`: bus active power injections.
- `bus_reactive_power_injections::Matrix{Float64}`: bus reactive power injections.
- `bus_active_power_withdrawals::Matrix{Float64}`: bus active power withdrawals.
- `bus_reactive_power_withdrawals::Matrix{Float64}`: bus reactive power withdrawals.
- `arc_active_power_flow_from_to::Matrix{Float64}`: active power flow measured at the from bus.
- `arc_reactive_power_flow_from_to::Matrix{Float64}`: reactive power flow measured at the from bus.
- `arc_active_power_flow_to_from::Matrix{Float64}`: active power flow measured at the to bus.
- `arc_reactive_power_flow_to_from::Matrix{Float64}`: reactive power flow measured at the to bus.
- `arc_angle_differences::Matrix{Float64}`: voltage angle difference across each arc.
- `converged::BitVector`: convergence status for each time step.
- `loss_factors::Union{Matrix{Float64}, Nothing}`: bus loss factors, or `nothing` if not computed.
- `voltage_stability_factors::Union{Matrix{Float64}, Nothing}`: voltage stability factors, or `nothing` if not computed.
- `arc_active_power_losses::Union{Matrix{Float64}, Nothing}`: active power losses per arc, or `nothing` if not computed.
- `lcc_active_power_flow_from_to::Union{Matrix{Float64}, Nothing}`: LCC active power flow at the from bus, or `nothing` if no LCCs.
- `lcc_active_power_flow_to_from::Union{Matrix{Float64}, Nothing}`: LCC active power flow at the to bus, or `nothing` if no LCCs.
- `lcc_reactive_power_flow_from_to::Union{Matrix{Float64}, Nothing}`: LCC reactive power flow at the from bus, or `nothing` if no LCCs.
- `lcc_reactive_power_flow_to_from::Union{Matrix{Float64}, Nothing}`: LCC reactive power flow at the to bus, or `nothing` if no LCCs.
"""
struct TimePowerFlowData <: AbstractPowerFlowResults
    bus_magnitude::Matrix{Float64}
    bus_angles::Matrix{Float64}
    bus_type::Matrix{PSY.ACBusTypes}
    bus_active_power_injections::Matrix{Float64}
    bus_reactive_power_injections::Matrix{Float64}
    bus_active_power_withdrawals::Matrix{Float64}
    bus_reactive_power_withdrawals::Matrix{Float64}
    arc_active_power_flow_from_to::Matrix{Float64}
    arc_reactive_power_flow_from_to::Matrix{Float64}
    arc_active_power_flow_to_from::Matrix{Float64}
    arc_reactive_power_flow_to_from::Matrix{Float64}
    arc_angle_differences::Matrix{Float64}
    converged::BitVector
    loss_factors::Union{Matrix{Float64}, Nothing}
    voltage_stability_factors::Union{Matrix{Float64}, Nothing}
    arc_active_power_losses::Union{Matrix{Float64}, Nothing}
    lcc_active_power_flow_from_to::Union{Matrix{Float64}, Nothing}
    lcc_active_power_flow_to_from::Union{Matrix{Float64}, Nothing}
    lcc_reactive_power_flow_from_to::Union{Matrix{Float64}, Nothing}
    lcc_reactive_power_flow_to_from::Union{Matrix{Float64}, Nothing}
end

"""
    TimePowerFlowData(n_buses::Int, n_arcs::Int, n_time_steps::Int;
        n_lccs::Int = 0,
        calculate_loss_factors::Bool = false,
        calculate_voltage_stability_factors::Bool = false,
        make_arc_active_power_losses::Bool = false)

Construct a `TimePowerFlowData` with pre-allocated arrays for the given dimensions.

Bus magnitudes default to 1.0 (flat start), bus angles default to 0.0, and bus types
default to `PSY.ACBusTypes.PQ`. Arc fields default to zeros. Optional fields
(`loss_factors`, `voltage_stability_factors`, `arc_active_power_losses`) are `nothing`
unless the corresponding keyword argument is `true`. LCC flow fields are allocated
when `n_lccs > 0`, otherwise `nothing`.
"""
function TimePowerFlowData(
    n_buses::Int,
    n_arcs::Int,
    n_time_steps::Int;
    n_lccs::Int = 0,
    calculate_loss_factors::Bool = false,
    calculate_voltage_stability_factors::Bool = false,
    make_arc_active_power_losses::Bool = false,
)
    return TimePowerFlowData(
        ones(n_buses, n_time_steps),
        zeros(n_buses, n_time_steps),
        fill(PSY.ACBusTypes.PQ, (n_buses, n_time_steps)),
        zeros(n_buses, n_time_steps),
        zeros(n_buses, n_time_steps),
        zeros(n_buses, n_time_steps),
        zeros(n_buses, n_time_steps),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        falses(n_time_steps),
        calculate_loss_factors ? zeros(n_buses, n_time_steps) : nothing,
        calculate_voltage_stability_factors ? zeros(n_buses, n_time_steps) : nothing,
        make_arc_active_power_losses ? zeros(n_arcs, n_time_steps) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps) : nothing,
    )
end

"""
    TimeContingencyPowerFlowData <: AbstractPowerFlowResults

Container for time-series contingency power flow results, storing bus and arc quantities
across multiple time steps and contingencies as 3D arrays with dimensions
`(entity, time_step, contingency)`.

# Fields
- `bus_magnitude::Array{Float64, 3}`: bus voltage magnitudes.
- `bus_angles::Array{Float64, 3}`: bus voltage angles.
- `bus_type::Array{PSY.ACBusTypes, 3}`: bus type classification at each time step and contingency.
- `bus_active_power_injections::Array{Float64, 3}`: bus active power injections.
- `bus_reactive_power_injections::Array{Float64, 3}`: bus reactive power injections.
- `bus_active_power_withdrawals::Array{Float64, 3}`: bus active power withdrawals.
- `bus_reactive_power_withdrawals::Array{Float64, 3}`: bus reactive power withdrawals.
- `arc_active_power_flow_from_to::Array{Float64, 3}`: active power flow measured at the from bus.
- `arc_reactive_power_flow_from_to::Array{Float64, 3}`: reactive power flow measured at the from bus.
- `arc_active_power_flow_to_from::Array{Float64, 3}`: active power flow measured at the to bus.
- `arc_reactive_power_flow_to_from::Array{Float64, 3}`: reactive power flow measured at the to bus.
- `arc_angle_differences::Array{Float64, 3}`: voltage angle difference across each arc.
- `converged::BitMatrix`: convergence status with dimensions `(time_steps, contingencies)`.
- `loss_factors::Union{Array{Float64, 3}, Nothing}`: bus loss factors, or `nothing` if not computed.
- `voltage_stability_factors::Union{Array{Float64, 3}, Nothing}`: voltage stability factors, or `nothing` if not computed.
- `arc_active_power_losses::Union{Array{Float64, 3}, Nothing}`: active power losses per arc, or `nothing` if not computed.
- `lcc_active_power_flow_from_to::Union{Array{Float64, 3}, Nothing}`: LCC active power flow at the from bus, or `nothing` if no LCCs.
- `lcc_active_power_flow_to_from::Union{Array{Float64, 3}, Nothing}`: LCC active power flow at the to bus, or `nothing` if no LCCs.
- `lcc_reactive_power_flow_from_to::Union{Array{Float64, 3}, Nothing}`: LCC reactive power flow at the from bus, or `nothing` if no LCCs.
- `lcc_reactive_power_flow_to_from::Union{Array{Float64, 3}, Nothing}`: LCC reactive power flow at the to bus, or `nothing` if no LCCs.
- `contingency_labels::Vector{String}`: human-readable labels for each contingency.
- `contingency_lookup::Dict{String, Int}`: mapping from contingency label to index.
- `network_modifications::Vector{Union{Nothing, PNM.NetworkModification}}`: network modification for each contingency.
"""
struct TimeContingencyPowerFlowData <: AbstractPowerFlowResults
    bus_magnitude::Array{Float64, 3}
    bus_angles::Array{Float64, 3}
    bus_type::Array{PSY.ACBusTypes, 3}
    bus_active_power_injections::Array{Float64, 3}
    bus_reactive_power_injections::Array{Float64, 3}
    bus_active_power_withdrawals::Array{Float64, 3}
    bus_reactive_power_withdrawals::Array{Float64, 3}
    arc_active_power_flow_from_to::Array{Float64, 3}
    arc_reactive_power_flow_from_to::Array{Float64, 3}
    arc_active_power_flow_to_from::Array{Float64, 3}
    arc_reactive_power_flow_to_from::Array{Float64, 3}
    arc_angle_differences::Array{Float64, 3}
    converged::BitMatrix
    loss_factors::Union{Array{Float64, 3}, Nothing}
    voltage_stability_factors::Union{Array{Float64, 3}, Nothing}
    arc_active_power_losses::Union{Array{Float64, 3}, Nothing}
    lcc_active_power_flow_from_to::Union{Array{Float64, 3}, Nothing}
    lcc_active_power_flow_to_from::Union{Array{Float64, 3}, Nothing}
    lcc_reactive_power_flow_from_to::Union{Array{Float64, 3}, Nothing}
    lcc_reactive_power_flow_to_from::Union{Array{Float64, 3}, Nothing}
    contingency_labels::Vector{String}
    contingency_lookup::Dict{String, Int}
    network_modifications::Vector{Union{Nothing, PNM.NetworkModification}}
    # Islanding reporting — populated by the contingency solve. A contingency
    # can island the network and still be "converged" in the sense that the
    # solve produced a minimum-norm pseudo-inverse solution (relative angles
    # within each island are meaningful). Absolute angles across islands are
    # not. Downstream consumers reading absolute angles MUST check `islanded`.
    islanded::Vector{Bool}
    n_islands::Vector{Int}
end

"""
    TimeContingencyPowerFlowData(n_buses::Int, n_arcs::Int, n_time_steps::Int,
        contingency_labels::Vector{String};
        n_lccs::Int = 0,
        network_modifications::Vector{Union{Nothing, PNM.NetworkModification}} = ...,
        calculate_loss_factors::Bool = false,
        calculate_voltage_stability_factors::Bool = false,
        make_arc_active_power_losses::Bool = false)

Construct a `TimeContingencyPowerFlowData` with pre-allocated 3D arrays of dimensions
`(entity, time_step, contingency)`.

Bus magnitudes default to 1.0 (flat start), bus angles default to 0.0, bus types default
to `PSY.ACBusTypes.PQ`, and arc fields default to zeros. Optional fields are `nothing`
unless the corresponding keyword argument is `true`. LCC flow fields are allocated
when `n_lccs > 0`, otherwise `nothing`. The `contingency_lookup` dictionary
is built automatically from `contingency_labels`.
"""
function TimeContingencyPowerFlowData(
    n_buses::Int,
    n_arcs::Int,
    n_time_steps::Int,
    contingency_labels::Vector{String};
    n_lccs::Int = 0,
    network_modifications::Vector{Union{Nothing, PNM.NetworkModification}} =
    Union{Nothing, PNM.NetworkModification}[nothing for _ in contingency_labels],
    calculate_loss_factors::Bool = false,
    calculate_voltage_stability_factors::Bool = false,
    make_arc_active_power_losses::Bool = false,
)
    n_ctg = length(contingency_labels)
    seen = Set{String}()
    duplicates = String[]
    for label in contingency_labels
        label in seen ? push!(duplicates, label) : push!(seen, label)
    end
    if !isempty(duplicates)
        throw(
            ArgumentError(
                "contingency_labels must be unique. Duplicates found: $(duplicates)",
            ),
        )
    end
    if length(network_modifications) != n_ctg
        throw(
            ArgumentError(
                "length(network_modifications) = $(length(network_modifications)) " *
                "must equal length(contingency_labels) = $n_ctg",
            ),
        )
    end
    contingency_lookup = Dict{String, Int}(
        label => idx for (idx, label) in enumerate(contingency_labels)
    )
    return TimeContingencyPowerFlowData(
        ones(n_buses, n_time_steps, n_ctg),
        zeros(n_buses, n_time_steps, n_ctg),
        fill(PSY.ACBusTypes.PQ, (n_buses, n_time_steps, n_ctg)),
        zeros(n_buses, n_time_steps, n_ctg),
        zeros(n_buses, n_time_steps, n_ctg),
        zeros(n_buses, n_time_steps, n_ctg),
        zeros(n_buses, n_time_steps, n_ctg),
        zeros(n_arcs, n_time_steps, n_ctg),
        zeros(n_arcs, n_time_steps, n_ctg),
        zeros(n_arcs, n_time_steps, n_ctg),
        zeros(n_arcs, n_time_steps, n_ctg),
        zeros(n_arcs, n_time_steps, n_ctg),
        falses(n_time_steps, n_ctg),
        calculate_loss_factors ? zeros(n_buses, n_time_steps, n_ctg) : nothing,
        calculate_voltage_stability_factors ?
        zeros(n_buses, n_time_steps, n_ctg) : nothing,
        make_arc_active_power_losses ? zeros(n_arcs, n_time_steps, n_ctg) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps, n_ctg) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps, n_ctg) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps, n_ctg) : nothing,
        n_lccs > 0 ? zeros(n_lccs, n_time_steps, n_ctg) : nothing,
        contingency_labels,
        contingency_lookup,
        network_modifications,
        fill(false, n_ctg),  # islanded — base (slot 1) is never islanded
        ones(Int, n_ctg),    # n_islands — default 1 (single connected component)
    )
end

# -- Accessor functions for TimeContingencyPowerFlowData --

"""Return the vector of contingency labels."""
get_contingency_labels(r::TimeContingencyPowerFlowData) = r.contingency_labels

"""Return the contingency label-to-index lookup dictionary."""
get_contingency_lookup(r::TimeContingencyPowerFlowData) = r.contingency_lookup

"""Return the vector of network modifications."""
get_network_modifications(r::TimeContingencyPowerFlowData) = r.network_modifications

"""Return the number of contingencies."""
get_n_contingencies(r::TimeContingencyPowerFlowData) = length(r.contingency_labels)

"""Return `true` if contingency `k` (or the one matching `label`) caused the
network to split into more than one electrical island. An islanded
contingency produces minimum-norm (pseudo-inverse) angles; ABSOLUTE angles
across distinct islands are not comparable, but relative angles within each
island are meaningful. `converged` remains `true` for islanded contingencies
that solved via pinv."""
get_islanded(r::TimeContingencyPowerFlowData, k::Int) = r.islanded[k]
get_islanded(r::TimeContingencyPowerFlowData, label::String) =
    r.islanded[r.contingency_lookup[label]]

"""Return the number of electrical islands after applying contingency `k` (or
`label`). `1` for a connected post-contingency network."""
get_n_islands(r::TimeContingencyPowerFlowData, k::Int) = r.n_islands[k]
get_n_islands(r::TimeContingencyPowerFlowData, label::String) =
    r.n_islands[r.contingency_lookup[label]]

"""
    get_network_modification(r::TimeContingencyPowerFlowData, label::String)

Return the `NetworkModification` (or `nothing`) for the contingency identified by `label`.
"""
function get_network_modification(r::TimeContingencyPowerFlowData, label::String)
    idx = r.contingency_lookup[label]
    return r.network_modifications[idx]
end

"""
    get_network_modification(r::TimeContingencyPowerFlowData, idx::Int)

Return the `NetworkModification` (or `nothing`) for the contingency at index `idx`.
"""
function get_network_modification(r::TimeContingencyPowerFlowData, idx::Int)
    return r.network_modifications[idx]
end

"""
    set_network_modification!(r::TimeContingencyPowerFlowData, label::String,
        mod::PNM.NetworkModification)

Set the network modification for the contingency identified by `label`.
"""
function set_network_modification!(
    r::TimeContingencyPowerFlowData,
    label::String,
    mod::PNM.NetworkModification,
)
    idx = r.contingency_lookup[label]
    r.network_modifications[idx] = mod
    return nothing
end

"""Result fields that are 3D arrays in `TimeContingencyPowerFlowData` and support slicing.
Derived from `_RESULT_FIELD_NAMES`, excluding `:converged` (which is a 2D `BitMatrix`)."""
const _CONTINGENCY_3D_FIELDS =
    filter(s -> s !== :converged, _RESULT_FIELD_NAMES)

"""
    get_contingency_slice(r::TimeContingencyPowerFlowData, field::Symbol, idx::Int)

Return a 2D view `(entity, time_step)` into the 3D result array for contingency `idx`.
Only 3D array fields are supported; optional fields that are `nothing` will throw.
"""
function get_contingency_slice(
    r::TimeContingencyPowerFlowData,
    field::Symbol,
    idx::Int,
)
    if !(field in _CONTINGENCY_3D_FIELDS)
        throw(
            ArgumentError(
                "field `$field` is not a valid 3D contingency result field. " *
                "Valid fields: $(sort(collect(_CONTINGENCY_3D_FIELDS)))",
            ),
        )
    end
    arr = getfield(r, field)
    if arr === nothing
        throw(
            ArgumentError(
                "field `$field` is `nothing` (was not computed). " *
                "Enable it at construction time.",
            ),
        )
    end
    n_ctg = get_n_contingencies(r)
    if idx < 1 || idx > n_ctg
        throw(
            ArgumentError(
                "contingency index $idx is out of bounds. Valid range: 1:$n_ctg",
            ),
        )
    end
    return @view arr[:, :, idx]
end

function get_contingency_slice(
    r::TimeContingencyPowerFlowData,
    field::Symbol,
    label::String,
)
    lookup = get_contingency_lookup(r)
    if !haskey(lookup, label)
        throw(
            ArgumentError(
                "contingency label \"$label\" not found. " *
                "Available labels: $(get_contingency_labels(r))",
            ),
        )
    end
    idx = lookup[label]
    return get_contingency_slice(r, field, idx)
end
