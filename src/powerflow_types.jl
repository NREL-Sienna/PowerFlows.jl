abstract type PowerFlowEvaluationModel end
abstract type ACPowerFlowSolverType end

struct KLUACPowerFlow <: ACPowerFlowSolverType end
struct NLSolveACPowerFlow <: ACPowerFlowSolverType end

Base.@kwdef struct ACPowerFlow{ACSolver <: ACPowerFlowSolverType} <:
                   PowerFlowEvaluationModel
    check_reactive_power_limits::Bool = false
end

# Create a constructor that defaults to KLUACPowerFlow
function ACPowerFlow(;
    check_reactive_power_limits::Bool = false,
    ACSolver::Type{<:ACPowerFlowSolverType} = KLUACPowerFlow,
)
    return ACPowerFlow{ACSolver}(check_reactive_power_limits)
end

struct DCPowerFlow <: PowerFlowEvaluationModel end
struct PTDFDCPowerFlow <: PowerFlowEvaluationModel end
struct vPTDFDCPowerFlow <: PowerFlowEvaluationModel end

Base.@kwdef struct PSSEExportPowerFlow <: PowerFlowEvaluationModel
    psse_version::Symbol
    export_dir::AbstractString
    write_comments::Bool = false
end
