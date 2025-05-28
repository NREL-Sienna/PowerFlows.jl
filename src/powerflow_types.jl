"""An abstract supertype for all types of power flows.
Subtypes: [`ACPowerFlow`](@ref) and [`AbstractDCPowerFlow`](@ref)."""
abstract type PowerFlowEvaluationModel end

"""An abstract supertype for all iterative methods.
Subtypes: [`NewtonRaphsonACPowerFlow`](@ref), [`TrustRegionACPowerFlow`](@ref), 
[`LevenbergMarquardtACPowerFlow`](@ref), and [`RobustHomotopyPowerFlow`](@ref).
"""
abstract type ACPowerFlowSolverType end

"""
    NewtonRaphsonACPowerFlow <: ACPowerFlowSolverType

An [`ACPowerFlowSolverType`](@ref) corresponding to a basic Newton-Raphson iterative method. 
The Newton step is taken verbatim at each iteration: no line search is performed.
"""
struct NewtonRaphsonACPowerFlow <: ACPowerFlowSolverType end

"""
    TrustRegionACPowerFlow <: ACPowerFlowSolverType

An [`ACPowerFlowSolverType`](@ref) corresponding to the Powell dogleg iterative method."""
struct TrustRegionACPowerFlow <: ACPowerFlowSolverType end

"""
    LevenbergMarquardtACPowerFlow <: ACPowerFlowSolverType

An [`ACPowerFlowSolverType`](@ref) corresponding to the Levenberg-Marquardt iterative method."""
struct LevenbergMarquardtACPowerFlow <: ACPowerFlowSolverType end

"""
    RobustHomotopyPowerFlow <: ACPowerFlowSolverType

An [`ACPowerFlowSolverType`](@ref) corresponding to a homotopy iterative method, based on the
paper [\"Improving the robustness of Newton-based power flow methods to cope with poor 
initial points\"](https://ieeexplore.ieee.org/document/6666905)."""
struct RobustHomotopyPowerFlow <: ACPowerFlowSolverType end

"""A struct for evaluating power flow solutions in AC systems.


This struct is parameterized by the type of AC power flow solver to use, which must be a
subtype of [`ACPowerFlowSolverType`](@ref). It also contains a few 
fields that control whether to compute certain additional data, like loss factors: 
see the constructor for details."""
struct ACPowerFlow{ACSolver <: ACPowerFlowSolverType} <: PowerFlowEvaluationModel
    check_reactive_power_limits::Bool
    exporter::Union{Nothing, PowerFlowEvaluationModel}
    calculate_loss_factors::Bool
    calculate_voltage_stability_factors::Bool
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
        Vector{Dict{Tuple{DataType, String}, Float64}},
    }
    enhanced_flat_start::Bool
    robust_power_flow::Bool
    skip_redistribution::Bool
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
    calculate_voltage_stability_factors::Bool = false,
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
        Vector{Dict{Tuple{DataType, String}, Float64}},
    } = nothing,
    enhanced_flat_start::Bool = true,
    robust_power_flow::Bool = false,
    skip_redistribution::Bool = false,
) where {ACSolver <: ACPowerFlowSolverType} =
    ACPowerFlow{ACSolver}(
        check_reactive_power_limits,
        exporter,
        calculate_loss_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        enhanced_flat_start,
        robust_power_flow,
        skip_redistribution,
    )

function ACPowerFlow(
    ACSolver::Type{<:ACPowerFlowSolverType} = NewtonRaphsonACPowerFlow;
    check_reactive_power_limits::Bool = false,
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    calculate_loss_factors::Bool = false,
    calculate_voltage_stability_factors::Bool = false,
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
        Vector{Dict{Tuple{DataType, String}, Float64}},
    } = nothing,
    enhanced_flat_start::Bool = true,
    robust_power_flow::Bool = false,
    skip_redistribution::Bool = false,
)
    if calculate_loss_factors && ACSolver == LevenbergMarquardtACPowerFlow
        error("Loss factor calculation is not supported by the Levenberg-Marquardt solver.")
    end

    return ACPowerFlow{ACSolver}(
        check_reactive_power_limits,
        exporter,
        calculate_loss_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        enhanced_flat_start,
        robust_power_flow,
        skip_redistribution,
    )
end

get_enhanced_flat_start(pf::ACPowerFlow) = pf.enhanced_flat_start
get_robust_power_flow(pf::ACPowerFlow) = pf.robust_power_flow
get_slack_participation_factors(pf::ACPowerFlow) = pf.generator_slack_participation_factors
get_calculate_loss_factors(pf::ACPowerFlow) = pf.calculate_loss_factors
get_calculate_voltage_stability_factors(pf::ACPowerFlow) =
    pf.calculate_voltage_stability_factors

abstract type AbstractDCPowerFlow <: PowerFlowEvaluationModel end

# only make sense for AC power flows, but convenient to have for code reuse reasons.
get_slack_participation_factors(::AbstractDCPowerFlow) = nothing
get_calculate_loss_factors(::AbstractDCPowerFlow) = false
get_calculate_voltage_stability_factors(::AbstractDCPowerFlow) = false

# the exporter field is not used in PowerFlows.jl, only in PowerSimulations.jl,
# which calls flatten_power_flow_evaluation_model then evaluates the two sequentially.
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
@kwdef struct DCPowerFlow <: AbstractDCPowerFlow
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
@kwdef struct PTDFDCPowerFlow <: AbstractDCPowerFlow
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
@kwdef struct vPTDFDCPowerFlow <: AbstractDCPowerFlow
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

"""
    PSSEExportPowerFlow(psse_version::Symbol, export_dir::AbstractString; kwargs...)

An evaluation model for exporting power flow results to PSSE format.

Arguments:
- `psse_version::Symbol`: The version of PSSE to export to. Must be among `$PSSE_EXPORT_SUPPORTED_VERSIONS`.
- `export_dir::AbstractString`: The directory where the PSSE files will be exported.
Optional keyword arguments:
- `name::AbstractString`: The base name for the exported files. Defaults to `\"$PSSE_DEFAULT_EXPORT_NAME\"`.
- `write_comments::Bool`: Whether to write comments in the exported files. Defaults to `false`.
- `overwrite::Bool`: Whether to overwrite the file if it exists already. Defaults to `false`.
"""
@kwdef struct PSSEExportPowerFlow <: PowerFlowEvaluationModel
    psse_version::Symbol
    export_dir::AbstractString
    name::AbstractString = PSSE_DEFAULT_EXPORT_NAME
    write_comments::Bool = false
    overwrite::Bool = false
end

PSSEExportPowerFlow(psse_version::Symbol, export_dir::AbstractString; kwargs...) =
    PSSEExportPowerFlow(; psse_version = psse_version, export_dir = export_dir, kwargs...)

get_exporter(pfem::PowerFlowEvaluationModel) = pfem.exporter
get_exporter(::PSSEExportPowerFlow) = nothing

"""
Expand a single `PowerFlowEvaluationModel` into its possibly multiple parts for separate
evaluation. Namely, if `pfem` contains a non-nothing `exporter`, return `[pfem, exporter]`,
else return `[pfem]`.
"""
function flatten_power_flow_evaluation_model(pfem::PowerFlowEvaluationModel)
    exporter = get_exporter(pfem)
    return if isnothing(exporter)
        PowerFlowEvaluationModel[pfem]
    else
        PowerFlowEvaluationModel[pfem, exporter]
    end
end
