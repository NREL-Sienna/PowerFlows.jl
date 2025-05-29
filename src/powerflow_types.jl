"""An abstract supertype for all types of power flows.

Subtypes: [`ACPowerFlow`](@ref), [`DCPowerFlow`](@ref), [`PTDFDCPowerFlow`](@ref), and [`vPTDFDCPowerFlow`](@ref)."""
abstract type PowerFlowEvaluationModel end

"""An abstract supertype for all iterative methods.

See [`NewtonRaphsonACPowerFlow`](@ref) and [`TrustRegionACPowerFlow`](@ref) for subtypes."""
abstract type ACPowerFlowSolverType end

"""
    NewtonRaphsonACPowerFlow <: ACPowerFlowSolverType

An [`ACPowerFlowSolverType`](@ref) corresponding to the Newton-Raphson iterative method."""
struct NewtonRaphsonACPowerFlow <: ACPowerFlowSolverType end
"""
    TrustRegionACPowerFlow <: ACPowerFlowSolverType

An [`ACPowerFlowSolverType`](@ref) corresponding to the Powell dogleg iterative method."""
struct TrustRegionACPowerFlow <: ACPowerFlowSolverType end

"""A struct for evaluating power flow solutions in AC systems.


This struct is parameterized by the type of AC power flow solver to use, which must be a
subtype of [`ACPowerFlowSolverType`](@ref). It also contains a few 
fields that control whether to compute certain additional data, like loss factors: 
see the constructor for details."""
struct ACPowerFlow{ACSolver <: ACPowerFlowSolverType} <: PowerFlowEvaluationModel
    check_reactive_power_limits::Bool
    exporter::Union{Nothing, PowerFlowEvaluationModel}
    calculate_loss_factors::Bool
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
        Vector{Dict{Tuple{DataType, String}, Float64}},
    }
end

"""
    ACPowerFlow{ACSolver}(
        check_reactive_power_limits::Bool = false,
        exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
        calculate_loss_factors::Bool = false,
        generator_slack_participation_factors::Union{
            Nothing,
            Dict{Tuple{DataType, String}, Float64},
            Vector{Dict{Tuple{DataType, String}, Float64}},
        } = nothing,
    ) where {ACSolver <: ACPowerFlowSolverType}

An evaluation model for a standard 
[AC powerflow](https://en.wikipedia.org/wiki/Power-flow_study#Power-flow_problem_formulation) 
with the specified solver type.


# Arguments
- `ACSolver`: The type of AC power flow solver to use, which must be a subtype of [`ACPowerFlowSolverType`](@ref).
    Default is [`NewtonRaphsonACPowerFlow`](@ref).
- `check_reactive_power_limits::Bool`: Whether to check reactive power limits during the power flow solution.
    Default is `false`.
- `exporter::Union{Nothing, PowerFlowEvaluationModel}`: An optional exporter for the power flow results. 
    If not `nothing`, it should be a [`PSSEExportPowerFlow`](@ref).
- `calculate_loss_factors::Bool`: Whether to calculate loss factors during the power flow solution.
    Default is `false`.
- `generator_slack_participation_factors::Union{Nothing, Dict{Tuple{DataType, String}, Float64}, Vector{Dict{Tuple{DataType, String}, Float64}}}`:
    An optional parameter that specifies the participation factors for generator slack in the power flow solution.
    If `nothing`, all slack is picked up by the reference bus. If a `Dict`, it should map `(component_type, component_name)`
    tuples to participation factors. If a `Vector`, it should contain multiple such dictionaries, 
    allowing for different participation factors for different time steps.
"""
ACPowerFlow{ACSolver}(;
    check_reactive_power_limits::Bool = false,
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    calculate_loss_factors::Bool = false,
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
        Vector{Dict{Tuple{DataType, String}, Float64}},
    } = nothing,
) where {ACSolver <: ACPowerFlowSolverType} =
    ACPowerFlow{ACSolver}(
        check_reactive_power_limits,
        exporter,
        calculate_loss_factors,
        generator_slack_participation_factors,
    )

ACPowerFlow(
    ACSolver::Type{<:ACPowerFlowSolverType} = NewtonRaphsonACPowerFlow;
    check_reactive_power_limits::Bool = false,
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    calculate_loss_factors::Bool = false,
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
        Vector{Dict{Tuple{DataType, String}, Float64}},
    } = nothing,
) = ACPowerFlow{ACSolver}(
    check_reactive_power_limits,
    exporter,
    calculate_loss_factors,
    generator_slack_participation_factors,
)

"""
    DCPowerFlow(
        exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    )

An evaluation model for a standard DC powerflow.

This provides a fast approximate solution to the AC powerflow problem, by solving for the 
bus voltage angles under some simplifying assumptions (lossless lines, constant voltage 
magnitudes, etc.). For details, see 
[Wikipedia](https://en.wikipedia.org/wiki/Power-flow_study#DC_power_flow)
or section 4 of the [MATPOWER docs](https://matpower.org/docs/MATPOWER-manual-4.1.pdf). If 
not `nothing`, the `exporter` should be a [`PSSEExportPowerFlow`](@ref).
"""
@kwdef struct DCPowerFlow <: PowerFlowEvaluationModel
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

"""
    PTDFDCPowerFlow(
        exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    )

An evaluation model that calculates line flows using the Power Transfer Distribution Factor 
Matrix.

This approximates the branch flows in the power grid, under some simplifying
assumptions (lossless lines, constant voltage magnitudes, etc.). See section 4 of the 
[MATPOWER docs](https://matpower.org/docs/MATPOWER-manual-4.1.pdf) for details. If not 
`nothing`, the `exporter` should be a [`PSSEExportPowerFlow`](@ref).
"""
@kwdef struct PTDFDCPowerFlow <: PowerFlowEvaluationModel
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

"""
    vPTDFDCPowerFlow(
        exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    )

An evaluation model that calculates line flows using a virtual Power Transfer Distribution 
Factor Matrix.

This is a replacement for the [PTDFDCPowerFlow](@ref) for large grids, 
where creating and storing the full PTDF matrix would be infeasible or slow. See the 
[PowerNetworkMatrices.jl docs](https://nrel-sienna.github.io/PowerNetworkMatrices.jl/stable/) for details. 
If not `nothing`, the `exporter` should be a [`PSSEExportPowerFlow`](@ref).
"""
@kwdef struct vPTDFDCPowerFlow <: PowerFlowEvaluationModel
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end
