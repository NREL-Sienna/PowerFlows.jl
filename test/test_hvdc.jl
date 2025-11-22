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

# written for below test case, but unused currently.
function add_component_with_power!(sys::PSY.System, bus::PSY.ACBus, P::Float64)
    if P > 0
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
    end
end

# goal: take a big system with a generic HVDC, replace its injections/withdrawals with
# generators, and verify that AC power flow converges to the same solution
# [except for the powers at the HVDC terminal buses]
# TODO currently errors, even on the unmodified system with the HVDCs.
# AC ignores the HVDC lines, so it sees multiple components, some without a REF bus.
#=
@testset "Test Generic HVDC on big network" begin
    sys_original = build_system(PSISystems, "HVDC_TWO_RTO_RTS_1Hr_sys")
    set_units_base_system!(sys_original, "SYSTEM_BASE")

    sys_modified = deepcopy(sys_original)
    set_units_base_system!(sys_modified, "SYSTEM_BASE")

    pf = ACPowerFlow{PF.TrustRegionACPowerFlow}()
    data_original = PF.PowerFlowData(
        pf,
        sys_original;
        correct_bustypes = true,
    )
    solve_powerflow!(data_original; pf = pf)
    println("original system solved successfully.")

    hvdc_bus_nums = Int[]
    for hvdc in get_components(PSY.TwoTerminalGenericHVDCLine, sys_modified)
        (P_from, P_to) = PF.hvdc_injections_natural_units(hvdc) .* PSY.get_base_power(sys_modified)
        P_from /= PSY.get_base_power(sys_modified)
        P_to /= PSY.get_base_power(sys_modified)
        arc = get_arc(hvdc)
        bus_from = arc.from
        bus_to = arc.to
        push!(hvdc_bus_nums, PSY.get_number(bus_from))
        push!(hvdc_bus_nums, PSY.get_number(bus_to))
        add_component_with_power!(sys_modified, bus_from, P_from)
        add_component_with_power!(sys_modified, bus_to, P_to)
        remove_component!(sys_modified, hvdc)
    end

    data_modified = PF.PowerFlowData(
        pf,
        sys_modified;
        correct_bustypes = true,
    )
    solve_powerflow!(data_modified; pf = pf)

    bus_lookup = PF.get_bus_lookup(data_original)
    hvdc_bus_inds = [bus_lookup[num] for num in hvdc_bus_nums]
    state_inds = [2 .* hvdc_bus_inds .- 1, 2 .* hvdc_bus_inds] # voltage mag and angle
    n_buses = length(PF.get_bus_axis(data_original))
    for i in 1:n_buses
        if i in hvdc_bus_inds
            continue
        end
        for fieldname in (:bus_magnitude, :bus_angle)
            @test isapprox(
                getfield(data_original, fieldname)[i],
                getfield(data_modified, fieldname)[i];
                atol = 1e-6,
            )
        end
        for fieldname in (:bus_activepower_injection, :bus_activepower_withdrawals,
                :bus_reactivepower_injection, :bus_reactivepower_withdrawals)
            @test isapprox(
                getfield(data_original, fieldname)[i],
                getfield(data_modified, fieldname)[i];
                atol = 1e-6,
            )
        end
    end
end
=#
