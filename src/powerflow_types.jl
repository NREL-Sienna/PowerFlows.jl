abstract type PowerFlowEvaluationModel end
abstract type ACPowerFlowSolverType end

struct KLUACPowerFlow <: ACPowerFlowSolverType end
struct NLSolveACPowerFlow <: ACPowerFlowSolverType end
struct LUACPowerFlow <: ACPowerFlowSolverType end  # Only for testing, a basic implementation using LinearAlgebra.lu, allocates a lot of memory

Base.@kwdef struct ACPowerFlow{ACSolver <: ACPowerFlowSolverType} <:
                   PowerFlowEvaluationModel
    check_reactive_power_limits::Bool = false
end

# Create a constructor for ACPowerFlow that defaults to KLUACPowerFlow
function ACPowerFlow(ACSolver::Type{<:ACPowerFlowSolverType} = KLUACPowerFlow;
    check_reactive_power_limits::Bool = false,
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
