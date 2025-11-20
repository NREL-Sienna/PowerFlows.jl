@testset "Test HVDC injections helper function" begin
    sys = build_system(MatpowerTestSystems, "matpower_case5_dc_sys")
    hvdc = only(get_components(TwoTerminalHVDC, sys))
    set_units_base_system!(sys, UnitSystem.NATURAL_UNITS)

    P_dc = 10.0
    set_active_power_flow!(hvdc, P_dc)

    # fixed slope loss function: 1% loss
    loss_coeff = 0.01
    set_loss!(hvdc, LinearCurve(loss_coeff))
    (P_from, P_to) = PF.hvdc_injections_natural_units(hvdc)
    @test isapprox(-P_from, P_dc)
    @test isapprox(P_to, P_dc * (1 - loss_coeff))

    # constant loss: 0.5 MW loss
    loss_const = 0.5
    set_loss!(hvdc, LinearCurve(0.0, loss_const))
    (P_from, P_to) = PF.hvdc_injections_natural_units(hvdc)
    @test isapprox(-P_from, P_dc)
    @test isapprox(P_to, P_dc - loss_const)

    # piecewise linear loss: 1% loss up to 20 MW, then 2% marginal loss above that
    loss_curve = PiecewiseIncrementalCurve(0.0, 0.0, [0.0, 20.0, 100.0], [0.01, 0.02])
    set_loss!(hvdc, loss_curve)
    for P_dc_setpoint in (10.0, 30.0)
        set_active_power_flow!(hvdc, P_dc_setpoint)
        (P_from, P_to) = PF.hvdc_injections_natural_units(hvdc)
        expected_loss =
            0.01 * min(P_dc_setpoint, 20.0) + 0.02 * max(0.0, P_dc_setpoint - 20.0)
        @test isapprox(-P_from, P_dc_setpoint)
        @test isapprox(P_to, P_dc_setpoint - expected_loss)
    end

    # test reversed flow error
    _, lcc = simple_lcc_system()
    set_active_power_flow!(lcc, -15.0)
    @test_throws ArgumentError PF.hvdc_injections_natural_units(lcc)
end

@testset "Test AC on generic HVDC" begin
    sys = build_system(MatpowerTestSystems, "matpower_case5_dc_sys")

    pf = ACPowerFlow{PF.TrustRegionACPowerFlow}()
    data = PF.PowerFlowData(
        pf,
        sys;
        correct_bustypes = true,
    )
    solve_powerflow!(data; pf = pf)
    @test all(data.converged)
end

@testset "Test DC on generic HVDC" begin
    for DC_type in (PF.DCPowerFlow, PF.PTDFDCPowerFlow, PF.vPTDFDCPowerFlow)
        sys = build_system(MatpowerTestSystems, "matpower_case5_dc_sys")

        pf = DC_type()
        data = PF.PowerFlowData(
            pf,
            sys;
        )
        solve_powerflow!(data)
        @test all(data.converged)
    end
end
