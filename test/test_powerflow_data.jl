@testset "PowerFlowData" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    PowerFlowData(ACPowerFlow(), sys)
    PowerFlowData(DCPowerFlow(), sys)
    PowerFlowData(PTDFDCPowerFlow(), sys)
    PowerFlowData(vPTDFDCPowerFlow(), sys)
end
