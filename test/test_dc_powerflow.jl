@testset "DC PowerFlow PTDF" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    data = PowerFlowData(PTDFDCPowerFlow(), sys)
    solve_powerflow!(data)
end

@testset "DC PowerFlow" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    data = PowerFlowData(DCPowerFlow(), sys)
    solve_powerflow!(data)
end
