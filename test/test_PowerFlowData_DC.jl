
@testset "test PowerFlowData structures" begin
    # load test sytem
    sys = PSB.build_system(PSITestSystems, "test_RTS_GMLC_sys")

    # get BA, ABA and PTDF matrix
    BA = PNM.BA_Matrix(sys)
    ABA = PNM.ABA_Matrix(sys; factorize=true)
    ptdf = PNM.PTDF(sys)

    # initialize PowerFlowData with PTDF matrix
    pfd_1 = PowerFlowData(DCPowerFlow(), sys, ptdf)

    # initialize PowerFlowData with ABA and BA matrices
    pfd_2 = PowerFlowData(DCPowerFlow(), sys, ABA, BA)

    # check flows between the different methods
    @test size(pfd_1.branch_flow_values) == size(pfd_2.branch_flow_values)

    for b in 1:pfd_1.n_branches
        @test pfd_1.branch_flow_values[b] == pfd_2.branch_flow_values[b]
    end
end