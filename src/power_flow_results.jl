"""
Abstract supertype for all power flow result containers.
"""
abstract type AbstractPowerFlowResults end

"""
    TimePowerFlowData <: AbstractPowerFlowResults

Container for time-series power flow results, storing bus and arc quantities across
multiple time steps as matrices (rows = buses/arcs, columns = time steps).

# Fields
- `bus_magnitude::Matrix{Float64}`: bus voltage magnitudes.
- `bus_angles::Matrix{Float64}`: bus voltage angles.
- `bus_type::Matrix{PSY.ACBusTypes}`: bus type classification at each time step.
- `arc_active_power_flow_from_to::Matrix{Float64}`: active power flow measured at the from bus.
- `arc_reactive_power_flow_from_to::Matrix{Float64}`: reactive power flow measured at the from bus.
- `arc_active_power_flow_to_from::Matrix{Float64}`: active power flow measured at the to bus.
- `arc_reactive_power_flow_to_from::Matrix{Float64}`: reactive power flow measured at the to bus.
- `arc_angle_differences::Matrix{Float64}`: voltage angle difference across each arc.
- `converged::BitVector`: convergence status for each time step.
- `loss_factors::Union{Matrix{Float64}, Nothing}`: bus loss factors, or `nothing` if not computed.
- `voltage_stability_factors::Union{Matrix{Float64}, Nothing}`: voltage stability factors, or `nothing` if not computed.
- `arc_active_power_losses::Union{Matrix{Float64}, Nothing}`: active power losses per arc, or `nothing` if not computed.
"""
struct TimePowerFlowData <: AbstractPowerFlowResults
    bus_magnitude::Matrix{Float64}
    bus_angles::Matrix{Float64}
    bus_type::Matrix{PSY.ACBusTypes}
    arc_active_power_flow_from_to::Matrix{Float64}
    arc_reactive_power_flow_from_to::Matrix{Float64}
    arc_active_power_flow_to_from::Matrix{Float64}
    arc_reactive_power_flow_to_from::Matrix{Float64}
    arc_angle_differences::Matrix{Float64}
    converged::BitVector
    loss_factors::Union{Matrix{Float64}, Nothing}
    voltage_stability_factors::Union{Matrix{Float64}, Nothing}
    arc_active_power_losses::Union{Matrix{Float64}, Nothing}
end

"""
    TimePowerFlowData(n_buses::Int, n_arcs::Int, n_time_steps::Int;
        calculate_loss_factors::Bool = false,
        calculate_voltage_stability_factors::Bool = false,
        make_arc_active_power_losses::Bool = false)

Construct a `TimePowerFlowData` with pre-allocated arrays for the given dimensions.

Bus magnitudes default to 1.0 (flat start), bus angles default to 0.0, and bus types
default to `PSY.ACBusTypes.PQ`. Arc fields default to zeros. Optional fields
(`loss_factors`, `voltage_stability_factors`, `arc_active_power_losses`) are `nothing`
unless the corresponding keyword argument is `true`.
"""
function TimePowerFlowData(
    n_buses::Int,
    n_arcs::Int,
    n_time_steps::Int;
    calculate_loss_factors::Bool = false,
    calculate_voltage_stability_factors::Bool = false,
    make_arc_active_power_losses::Bool = false,
)
    return TimePowerFlowData(
        ones(n_buses, n_time_steps),
        zeros(n_buses, n_time_steps),
        fill(PSY.ACBusTypes.PQ, (n_buses, n_time_steps)),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        zeros(n_arcs, n_time_steps),
        falses(n_time_steps),
        calculate_loss_factors ? zeros(n_buses, n_time_steps) : nothing,
        calculate_voltage_stability_factors ? zeros(n_buses, n_time_steps) : nothing,
        make_arc_active_power_losses ? zeros(n_arcs, n_time_steps) : nothing,
    )
end
