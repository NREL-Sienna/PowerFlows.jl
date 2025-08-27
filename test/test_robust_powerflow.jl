@testset "test robust homotopy powerflow" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    sys2 = deepcopy(sys)
    pf_hom = ACPowerFlow{PF.RobustHomotopyPowerFlow}()
    data_hom = PowerFlowData(pf_hom, sys)
    # infologger = ConsoleLogger(stderr, Logging.Info)
    # with_logger(infologger) do; solve_powerflow!(data_hom; pf = pf_hom); end;
    solve_powerflow!(data_hom; pf = pf_hom)

    pf_nr = ACPowerFlow()
    data_nr = PowerFlowData(pf_nr, sys2)
    solve_powerflow!(data_nr; pf = pf_nr)
    @test data_nr.bus_angles ≈ data_hom.bus_angles
    @test data_nr.bus_magnitude ≈ data_hom.bus_magnitude
end
