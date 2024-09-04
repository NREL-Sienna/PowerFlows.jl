abstract type PowerFlowEvaluationModel end

Base.@kwdef struct ACPowerFlow <: PowerFlowEvaluationModel
    check_reactive_power_limits::Bool = false
end

struct DCPowerFlow <: PowerFlowEvaluationModel end
struct PTDFDCPowerFlow <: PowerFlowEvaluationModel end
struct vPTDFDCPowerFlow <: PowerFlowEvaluationModel end

Base.@kwdef struct PSSEExportPowerFlow <: PowerFlowEvaluationModel
    psse_version::Symbol
    export_dir::AbstractString
    write_comments::Bool = false
end
