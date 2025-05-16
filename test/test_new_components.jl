# FIXME this test fails. The interface for some of these objects has changed.
@testset "new psy 5 components" begin
    # TODO are there other components that should be included here?
    componentToSys = Dict{Type, Tuple{Type, String}}(
        TapTransformer => (MatpowerTestSystems, "matpower_case2_sys"),
        TModelHVDCLine => (PSISystems, "sys10_pjm_ac_dc"),
        TwoTerminalVSCLine => (PSSEParsingTestSystems, "pti_vsc_hvdc_test_sys"),
        SwitchedAdmittance => (PSSEParsingTestSystems, "pti_case24_sys"),
        TwoTerminalLCCLine => (PSSEParsingTestSystems, "pti_two_terminal_hvdc_test_sys"),
        DiscreteControlledACBranch => (PSSEParsingTestSystems, "pti_modified_case14_sys"),
        TwoTerminalGenericHVDCLine => (PSITestSystems, "c_sys5_dc"),
        Transformer2W => (PSIDSystems, "psid_4bus_multigen"),
        Transformer3W => (PSSEParsingTestSystems, "pti_frankenstein_20_sys"),
        FACTSControlDevice => (PSSEParsingTestSystems, "pti_frankenstein_70_sys"),
        PhaseShiftingTransformer => (PSSEParsingTestSystems, "pti_two_winding_mag_test_sys")
    )
    for (component, sysInfo) in componentToSys
        sys = build_system(sysInfo...)
        pf = ACPowerFlow()
        data = PowerFlowData(pf, sys)
        @test solve_powerflow(pf, sys)
    end
end