abstract type PowerFlowEvaluationModel end
abstract type ACPowerFlowSolverType end

struct NewtonRaphsonACPowerFlow <: ACPowerFlowSolverType end
struct TrustRegionACPowerFlow <: ACPowerFlowSolverType end

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
    robust_power_flow::Bool
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
    robust_power_flow::Bool = false,
) where {ACSolver <: ACPowerFlowSolverType} =
    ACPowerFlow{ACSolver}(
        check_reactive_power_limits,
        exporter,
        calculate_loss_factors,
        calculate_voltage_stability_factors,
        generator_slack_participation_factors,
        robust_power_flow,
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
    robust_power_flow::Bool = false,
) = ACPowerFlow{ACSolver}(
    check_reactive_power_limits,
    exporter,
    calculate_loss_factors,
    calculate_voltage_stability_factors,
    generator_slack_participation_factors,
    robust_power_flow,
)

get_robust_power_flow(pf::ACPowerFlow{ACSolver}) where {ACSolver} = pf.robust_power_flow

abstract type AbstractDCPowerFlow <: PowerFlowEvaluationModel end

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
