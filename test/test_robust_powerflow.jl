@testset "test robust homotopy powerflow" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    pf_hom = ACPowerFlow{RobustHomoptyPowerFlow}()
    data_hom = PowerFlowData(pf_hom, sys)
    solve_powerflow!(data_hom; pf = pf_hom)
end