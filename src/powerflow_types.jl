Base.@kwdef struct ACPowerFlow
    check_reactive_power_limits::Bool = false
end

struct DCPowerFlow end
struct PTDFDCPowerFlow end
struct vPTDFDCPowerFlow end
