abstract type PowerFlowEvaluationModel end
abstract type ACPowerFlowSolverType end

struct KLUACPowerFlow <: ACPowerFlowSolverType end
struct HybridACPowerFlow <: ACPowerFlowSolverType end

struct ACPowerFlow{ACSolver <: ACPowerFlowSolverType} <: PowerFlowEvaluationModel
    check_reactive_power_limits::Bool
    exporter::Union{Nothing, PowerFlowEvaluationModel}
    calc_loss_factors::Bool
end

ACPowerFlow{ACSolver}(;
    check_reactive_power_limits::Bool = false,
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    calc_loss_factors::Bool = false,
) where {ACSolver <: ACPowerFlowSolverType} =
    ACPowerFlow{ACSolver}(check_reactive_power_limits, exporter, calc_loss_factors)

ACPowerFlow(
    ACSolver::Type{<:ACPowerFlowSolverType} = KLUACPowerFlow;
    check_reactive_power_limits::Bool = false,
    exporter::Union{Nothing, PowerFlowEvaluationModel} = nothing,
    calc_loss_factors::Bool = false,
) = ACPowerFlow{ACSolver}(check_reactive_power_limits, exporter, calc_loss_factors)

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
