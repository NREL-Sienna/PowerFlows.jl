@testset "test robust homotopy powerflow" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    pf_hom = ACPowerFlow{PF.RobustHomotopyPowerFlow}()
    data_hom = PowerFlowData(pf_hom, sys)
    time_step = 1
    x0 = PF.calculate_x0(data_hom, time_step)
    h = PF.HomotopyHessian(data_hom, time_step)

    solve_powerflow!(data_hom; pf = pf_hom)
end