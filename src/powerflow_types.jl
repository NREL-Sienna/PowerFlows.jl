abstract type PowerFlowEvaluationModel end
abstract type ACPowerFlowSolverType end

# note: refinement strategy and parameters are very basic, not tuned at all.
# so while those parameters are exposed, changing probably won't do much.
@kwdef struct NewtonRaphsonACPowerFlow <: ACPowerFlowSolverType
    max_iterations::Int = DEFAULT_NR_MAX_ITER
    tolerance::Float64 = DEFAULT_NR_TOL
    refinement_threshold::Float64 = DEFAULT_REFINEMENT_THRESHOLD
    refinement_epsilon::Float64 = DEFAULT_REFINEMENT_EPS
end

@kwdef struct TrustRegionACPowerFlow <: ACPowerFlowSolverType
    max_iterations::Int = DEFAULT_NR_MAX_ITER
    tolerance::Float64 = DEFAULT_NR_TOL
    refinement_threshold::Float64 = DEFAULT_REFINEMENT_THRESHOLD
    refinement_epsilon::Float64 = DEFAULT_REFINEMENT_EPS
    eta::Float64 = DEFAULT_TRUST_REGION_ETA
    trust_region_factor::Float64 = DEFAULT_TRUST_REGION_FACTOR
end

struct ACPowerFlow{ACSolver <: ACPowerFlowSolverType} <: PowerFlowEvaluationModel
    solver::ACSolver
    check_reactive_power_limits::Bool
    exporter::Union{Nothing, PowerFlowEvaluationModel}
    calculate_loss_factors::Bool
    generator_slack_participation_factors::Union{
        Nothing,
        Dict{Tuple{DataType, String}, Float64},
        Vector{Dict{Tuple{DataType, String}, Float64}},
    }
end

# fully specified solver instance.
ACPowerFlow(
    solver::ACSolver;
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
        solver,
        check_reactive_power_limits,
        exporter,
        calculate_loss_factors,
        generator_slack_participation_factors,
    )

# only type of solver provided: construct instance with default parameters. 
ACPowerFlow{ACSolver}(; kwargs...) where {ACSolver <: ACPowerFlowSolverType} =
    ACPowerFlow(ACSolver(); kwargs...)

# nothing specified: construct with default NR solver and default parameters.
ACPowerFlow(; kwargs...) = ACPowerFlow{NewtonRaphsonACPowerFlow}(; kwargs...)

@kwdef struct DCPowerFlow <: PowerFlowEvaluationModel
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

@kwdef struct PTDFDCPowerFlow <: PowerFlowEvaluationModel
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

@kwdef struct vPTDFDCPowerFlow <: PowerFlowEvaluationModel
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

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

get_solver(pfem::ACPowerFlow) = pfem.solver

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
