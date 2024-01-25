@testset "PowerFlowData" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    PowerFlowData(ACPowerFlow(), sys)
    PowerFlowData(DCPowerFlow(), sys)
    PowerFlowData(PTDFDCPowerFlow(), sys)
    PowerFlowData(vPTDFDCPowerFlow(), sys)
end

@testset "PowerFlowData multiperiod" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    time_steps = 24
    # TODO: "multiperiod AC still to implement"
    PowerFlowData(DCPowerFlow(), sys; time_steps = time_steps)
    PowerFlowData(PTDFDCPowerFlow(), sys; time_steps = time_steps)
    PowerFlowData(vPTDFDCPowerFlow(), sys; time_steps = time_steps)
end
