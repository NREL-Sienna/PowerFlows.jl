@testset "PowerFlowData" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    PowerFlowData(ACPowerFlow(), sys)
    PowerFlowData(DCPowerFlow(), sys)
    PowerFlowData(PTDFDCPowerFlow(), sys)
    PowerFlowData(vPTDFDCPowerFlow(), sys)
end

@testset "PowerFlowData multiperiod" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    timesteps = 24
    # TODO: "multiperiod AC still to implement"
    PowerFlowData(DCPowerFlow(), sys; timesteps=timesteps)
    PowerFlowData(PTDFDCPowerFlow(), sys; timesteps=timesteps)
    PowerFlowData(vPTDFDCPowerFlow(), sys; timesteps=timesteps)
end
