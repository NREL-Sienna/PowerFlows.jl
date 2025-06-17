@testset "flow_func" begin
    b = PSY.PhaseShiftingTransformer(nothing)
    @test_throws ErrorException PF.flow_func(b, zero(ComplexF64), zero(ComplexF64))
end

# known issue.
@testset "_get_load_data" begin
    sys = System(100.0)
    b1 = ACBus(;
        number = 9,
        name = "bus_1",
        available = true,
        bustype = PSY.ACBusTypes.PQ,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b1)
    g1_active, g1_reactive = 1.0, 0.1
    g1 = StandardLoad(;
        name = "load_1",
        available = true,
        bus = b1,
        base_power = 100.0,
        constant_active_power = g1_active,
        constant_reactive_power = g1_reactive,
    )
    g2_active, g2_reactive = 0.8, 0.2
    g2 = StandardLoad(;
        name = "load_2",
        available = true,
        bus = b1,
        base_power = 100.0,
        constant_active_power = g2_active,
        constant_reactive_power = g2_reactive,
    )
    add_component!(sys, g1)
    add_component!(sys, g2)
    for (available1, available2) in Iterators.product([true, false], [true, false])
        set_available!(g1, available1)
        @assert g1.available == available1
        set_available!(g2, available2)
        @assert g2.available == available2
        active_total, reactive_total = PF._get_load_data(sys, b1)
        active_total_true = Int(available1) * g1_active + Int(available2) * g2_active
        reactive_total_true = Int(available1) * g1_reactive + Int(available2) * g2_reactive
        @test active_total == active_total_true
        @test reactive_total == reactive_total_true
    end
end

@testset "_power_redistribution_ref" begin
    
    # edge case: no devices at the ref bus.
    sys = System(100.0)
    b1 = ACBus(;
        number = 9,
        name = "bus_1",
        available = true,
        bustype = PSY.ACBusTypes.REF,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b1)
    @test_throws ErrorException PF._power_redistribution_ref(sys, 0.0, 0.0, b1, 1)

    # TODO get a proper test case here
    # 95% of this is written by ChatGPT: it's catered to the code itself,
    # not based on its expected behavior.
    for genA_activepower_limits in
                              [(min = 0.0, max = 12.0), (min = 0.0, max = 10.0)]
        sys = System(100.0)


        bus = ACBus(;
            number = 1,
            name = "Bus1",
            available = true,
            bustype = PSY.ACBusTypes.PV,
            angle = 0.0,
            magnitude = 1.0,
            voltage_limits = (0.0, 2.0),
            base_voltage = 230.0,
        )
        add_component!(sys, bus)

        # Define thermal generators with constrained P limits
        gen_a = ThermalStandard(;
            name = "GenA",
            available = true,
            status = true,
            bus = bus,
            active_power = 0.0,
            reactive_power = 0.0,
            rating = 10.0,
            base_power = 100.0,
            prime_mover_type = PrimeMovers.ST,
            fuel = ThermalFuels.COAL,
            active_power_limits = genA_activepower_limits,
            reactive_power_limits = (min = -5.0, max = 5.0),
            ramp_limits = (up = 0.2, down = 0.2),
            operation_cost = ThermalGenerationCost(nothing),
            time_limits = (up = 8.0, down = 8.0),
            must_run = false,
        )

        gen_b = ThermalStandard(;
            name = "GenB",
            available = true,
            status = true,
            bus = bus,
            active_power = 0.0,
            reactive_power = 0.0,
            rating = 5.0,
            base_power = 100.0,
            prime_mover_type = PrimeMovers.ST,
            fuel = ThermalFuels.NATURAL_GAS,
            active_power_limits = (min = 0.0, max = 5.0),
            reactive_power_limits = (min = -5.0, max = 5.0),
            ramp_limits = (up = 0.2, down = 0.2),
            operation_cost = ThermalGenerationCost(nothing),
            time_limits = (up = 8.0, down = 8.0),
            must_run = false,
        )

        gen_c = ThermalStandard(;
            name = "GenC",
            available = true,
            status = true,
            bus = bus,
            active_power = 0.0,
            reactive_power = 0.0,
            rating = 3.0,
            base_power = 100.0,
            prime_mover_type = PrimeMovers.ST,
            fuel = ThermalFuels.NUCLEAR,
            active_power_limits = (min = 0.0, max = 3.0),
            reactive_power_limits = (min = -5.0, max = 5.0),
            ramp_limits = (up = 0.2, down = 0.2),
            operation_cost = ThermalGenerationCost(nothing),
            time_limits = (up = 8.0, down = 8.0),
            must_run = false,
        )

        # Add components to system
        add_component!(sys, gen_a)
        add_component!(sys, gen_b)
        add_component!(sys, gen_c)

        PF._power_redistribution_ref(sys, 20.0, 0.0, bus, 5)
    end
end
