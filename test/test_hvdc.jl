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

function add_component_with_power!(sys::PSY.System, bus::PSY.ACBus, P::Float64)
    if P == 0.0
        println(
            "Skipping addition of component with zero power at bus $(PSY.get_number(bus))",
        )
        return
    elseif P > 0
        gen = ThermalStandard(;
            name = "gen_$(PSY.get_number(bus))_thermal_hvdc_$P",
            available = true,
            status = true,
            bus = bus,
            active_power = P, # TODO check units
            reactive_power = 0.0,
            # rest shouldn't actually matter for our test.
            rating = 1.0,
            active_power_limits = (min = 0, max = P),
            reactive_power_limits = (min = -1, max = 1),
            ramp_limits = nothing,
            operation_cost = ThermalGenerationCost(nothing),
            base_power = 100.0,
            time_limits = nothing,
            prime_mover_type = PrimeMovers.OT,
            fuel = ThermalFuels.OTHER,
            services = Device[],
            dynamic_injector = nothing,
            ext = Dict{String, Any}(),
        )
        add_component!(sys, gen)
        @assert get_active_power(gen) == P
        @assert PSY.get_units_setting(gen).unit_system == PSY.UnitSystem.SYSTEM_BASE
    else
        load = PowerLoad(;
            name = "load_$(PSY.get_number(bus))_hvdc_$(-P)",
            available = true,
            bus = bus,
            active_power = -P, # Per-unitized by device base_power
            reactive_power = 0.0, # Per-unitized by device base_power
            base_power = 100.0, # MVA
            max_active_power = -P,
            max_reactive_power = 0.0,
        )
        add_component!(sys, load)
        @assert get_active_power(load) == -P
        @assert PSY.get_units_setting(load).unit_system == PSY.UnitSystem.SYSTEM_BASE
    end
end

function replace_generic_hvdcs!(sys)
    for hvdc in get_components(PSY.TwoTerminalGenericHVDCLine, sys)
        if PSY.get_active_power_flow(hvdc) == 0.0
            @assert PSY.get_units_setting(hvdc).unit_system ==
                    PSY.UnitSystem.SYSTEM_BASE
            set_active_power_flow!(hvdc, 0.1)
        end
        (P_from, P_to) = PF.hvdc_injections_natural_units(hvdc)
        P_from /= PSY.get_base_power(sys)
        P_to /= PSY.get_base_power(sys)
        arc = get_arc(hvdc)
        bus_from = arc.from
        bus_to = arc.to
        add_component_with_power!(sys, bus_from, P_from)
        add_component_with_power!(sys, bus_to, P_to)
        remove_component!(sys, hvdc)
    end
    return sys
end

function test_generic_hvdc_on_big_system(pf_type::Type{<:PF.PowerFlowEvaluationModel})
    sys_original = build_system(PSISystems, "HVDC_TWO_RTO_RTS_1Hr_sys")
    set_units_base_system!(sys_original, "SYSTEM_BASE")

    for hvdc in get_components(PSY.TwoTerminalGenericHVDCLine, sys_original)
        @assert PSY.get_units_setting(hvdc).unit_system ==
                PSY.UnitSystem.SYSTEM_BASE
        set_active_power_flow!(hvdc, 0.1)
    end

    pf = pf_type()
    data_original = PF.PowerFlowData(
        pf,
        sys_original;
        correct_bustypes = true,
    )
    solve_powerflow!(data_original)
    ref_bus_inds = findall(data_original.bus_type .== (PSY.ACBusTypes.REF,))
    @test all(data_original.bus_angles[ref_bus_inds] .== 0.0)

    sys_modified = deepcopy(sys_original)
    set_units_base_system!(sys_modified, "SYSTEM_BASE")
    replace_generic_hvdcs!(sys_modified)

    data_modified = PF.PowerFlowData(
        pf,
        sys_modified;
        correct_bustypes = true,
    )
    solve_powerflow!(data_modified)

    # verify assumptions
    @assert !isempty(get_components(PSY.TwoTerminalGenericHVDCLine, sys_original))
    @assert isempty(get_components(PSY.TwoTerminalGenericHVDCLine, sys_modified))
    # arc axes aren't the same: the HVDC arcs must somehow change control flow in PNM.
    # @assert all(PF.get_arc_axis(data_original) .== PF.get_arc_axis(data_modified))
    @assert all(PF.get_bus_axis(data_original) .== PF.get_bus_axis(data_modified))

    @test all(
        abs.(
            data_original.power_network_matrix.data .-
            data_modified.power_network_matrix.data
        ) .< TIGHT_TOLERANCE,
    )
    if pf_type() isa PF.AbstractDCPowerFlow
        # for AC, aux network matrix defaults to nothing.
        @test all(
            abs.(
                data_original.aux_network_matrix.data .-
                data_modified.aux_network_matrix.data
            ) .< TIGHT_TOLERANCE,
        )
        # no need to check arc flows: since angles and aux matrices match, flows match.
        # check that the b in our Ax = b is the same
        rhs_original = (
            data_original.bus_activepower_injection .-
            data_original.bus_activepower_withdrawals .+
            data_original.bus_hvdc_net_power
        )
        rhs_modified = (
            data_modified.bus_activepower_injection .-
            data_modified.bus_activepower_withdrawals
        )
        @test all(abs.(rhs_original - rhs_modified) .< TIGHT_TOLERANCE)
    end

    for fieldname in (:bus_angles, :bus_magnitude)
        @test all(
            abs.(
                getfield(data_original, fieldname) .-
                getfield(data_modified, fieldname)
            ) .< TIGHT_TOLERANCE,
        )
    end

    for fieldname in (:arc_activepower_flow_from_to, :arc_activepower_flow_to_from,
        :arc_reactivepower_flow_from_to, :arc_reactivepower_flow_to_from)
        for (original_arc_ind, arc) in enumerate(PF.get_arc_axis(data_original))
            modified_arc_ind = PF.get_arc_lookup(data_modified)[arc]
            @test isapprox(
                getfield(data_original, fieldname)[original_arc_ind],
                getfield(data_modified, fieldname)[modified_arc_ind],
                atol = TIGHT_TOLERANCE,
            )
        end
    end
end

# take a big system with a generic HVDC, replace its injections/withdrawals with
# generators, and verify that a power flow converges to the same solution
@testset "Test Generic HVDC on big network: AC power flow" begin
    test_generic_hvdc_on_big_system(PF.ACPowerFlow)
end

@testset "Test Generic HVDC on big network: DC power flow" begin
    test_generic_hvdc_on_big_system(PF.DCPowerFlow)
end
