@testset "DC PowerFlow PTDF" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    solve_powerflow!(PTDFDCPowerFlow(), sys)
end

@testset "DC PowerFlow" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    solve_powerflow!(DCPowerFlow(), sys)
end
