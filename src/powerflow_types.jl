abstract type PowerFlowEvaluationModel end
abstract type ACPowerFlowSolverType end

struct NewtonRaphsonACPowerFlow <: ACPowerFlowSolverType end
struct TrustRegionACPowerFlow <: ACPowerFlowSolverType end
struct LevenbergMarquardtACPowerFlow <: ACPowerFlowSolverType end
struct RobustHomotopyPowerFlow <: ACPowerFlowSolverType end

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

ACPowerFlow(
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
) = ACPowerFlow{ACSolver}(
    check_reactive_power_limits,
    exporter,
    calculate_loss_factors,
    calculate_voltage_stability_factors,
    generator_slack_participation_factors,
    enhanced_flat_start,
    robust_power_flow,
    skip_redistribution,
)

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
@kwdef struct DCPowerFlow <: AbstractDCPowerFlow
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

@kwdef struct PTDFDCPowerFlow <: AbstractDCPowerFlow
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing
end

@kwdef struct vPTDFDCPowerFlow <: AbstractDCPowerFlow
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
