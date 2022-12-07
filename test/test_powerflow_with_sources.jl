@testset "Multiple sources at ref" begin
    sys = System(100.0)
    b = Bus(
        number=1,
        name="01",
        bustype=BusTypes.REF,
        angle=0.0,
        magnitude=1.1,
        voltage_limits=(0.0, 2.0),
        base_voltage=230,
    )
    add_component!(sys, b)

    #Test two sources with equal and opposite P and Q
    s1 = Source(
        name="source_1",
        available=true,
        bus=b,
        active_power=0.5,
        reactive_power=0.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s1)
    s2 = Source(
        name="source_2",
        available=true,
        bus=b,
        active_power=-0.5,
        reactive_power=-0.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s2)
    @test run_powerflow!(sys, finite_diff=true)

    #Create power mismatch, test for error 
    set_active_power!(get_component(Source, sys, "source_1"), -0.4)
    @test_throws ErrorException(
        "Sources do not match P and/or Q requirements for reference bus.",
    ) run_powerflow!(sys, finite_diff=true)
end

@testset "Multiple sources at PV" begin
    sys = System(100.0)
    b1 = Bus(
        number=1,
        name="01",
        bustype=BusTypes.REF,
        angle=0.0,
        magnitude=1.1,
        voltage_limits=(0.0, 2.0),
        base_voltage=230,
    )
    add_component!(sys, b1)
    b2 = Bus(
        number=2,
        name="02",
        bustype=BusTypes.PV,
        angle=0.0,
        magnitude=1.1,
        voltage_limits=(0.0, 2.0),
        base_voltage=230,
    )
    add_component!(sys, b2)
    a = Arc(from=b1, to=b2)
    add_component!(sys, a)
    l = Line(
        name="l1",
        available=true,
        active_power_flow=0.0,
        reactive_power_flow=0.0,
        arc=a,
        r=1e-3,
        x=1e-3,
        b=(from=0.0, to=0.0),
        rate=0.0,
        angle_limits=(min=-pi / 2, max=pi / 2),
    )
    add_component!(sys, l)

    #Test two sources with equal and opposite P and Q
    s1 = Source(
        name="source_1",
        available=true,
        bus=b1,
        active_power=0.5,
        reactive_power=0.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s1)
    s2 = Source(
        name="source_2",
        available=true,
        bus=b2,
        active_power=0.5,
        reactive_power=1.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s2)
    s3 = Source(
        name="source_3",
        available=true,
        bus=b2,
        active_power=-0.5,
        reactive_power=-1.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s3)

    @test run_powerflow!(sys, finite_diff=true)

    #Create power mismatch, test for error 
    set_reactive_power!(get_component(Source, sys, "source_3"), -0.5)
    @test_throws ErrorException("Sources do not match Q requirements for PV bus.") run_powerflow!(
        sys,
        finite_diff=true,
    )
end

@testset "Source + non-source at Ref" begin
    sys = System(100.0)
    b = Bus(
        number=1,
        name="01",
        bustype=BusTypes.REF,
        angle=0.0,
        magnitude=1.1,
        voltage_limits=(0.0, 2.0),
        base_voltage=230,
    )
    add_component!(sys, b)

    #Test two sources with equal and opposite P and Q
    s1 = Source(
        name="source_1",
        available=true,
        bus=b,
        active_power=0.5,
        reactive_power=0.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s1)
    g1 = ThermalStandard(
        name="init",
        available=true,
        status=false,
        bus=b,
        active_power=0.1,
        reactive_power=0.1,
        rating=0.0,
        active_power_limits=(min=0.0, max=0.0),
        reactive_power_limits=nothing,
        ramp_limits=nothing,
        operation_cost=ThreePartCost(nothing),
        base_power=100.0,
        time_limits=nothing,
        prime_mover=PrimeMovers.OT,
        fuel=ThermalFuels.OTHER,
        services=Device[],
        dynamic_injector=nothing,
        ext=Dict{String, Any}(),
        time_series_container=InfrastructureSystems.TimeSeriesContainer(),
    )
    add_component!(sys, g1)

    @test run_powerflow!(sys, finite_diff=true)
    @test isapprox(
        get_active_power(get_component(Source, sys, "source_1")),
        0.5;
        atol=0.001,
    )
    @test isapprox(
        get_reactive_power(get_component(Source, sys, "source_1")),
        0.1;
        atol=0.001,
    )
end

@testset "Source + non-source at PV" begin
    sys = System(100.0)
    b1 = Bus(
        number=1,
        name="01",
        bustype=BusTypes.REF,
        angle=0.0,
        magnitude=1.1,
        voltage_limits=(0.0, 2.0),
        base_voltage=230,
    )
    add_component!(sys, b1)
    b2 = Bus(
        number=2,
        name="02",
        bustype=BusTypes.PV,
        angle=0.0,
        magnitude=1.1,
        voltage_limits=(0.0, 2.0),
        base_voltage=230,
    )
    add_component!(sys, b2)
    a = Arc(from=b1, to=b2)
    add_component!(sys, a)
    l = Line(
        name="l1",
        available=true,
        active_power_flow=0.0,
        reactive_power_flow=0.0,
        arc=a,
        r=1e-3,
        x=1e-3,
        b=(from=0.0, to=0.0),
        rate=0.0,
        angle_limits=(min=-pi / 2, max=pi / 2),
    )
    add_component!(sys, l)

    #Test two sources with equal and opposite P and Q
    s1 = Source(
        name="source_1",
        available=true,
        bus=b1,
        active_power=0.5,
        reactive_power=0.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s1)
    s2 = Source(
        name="source_2",
        available=true,
        bus=b2,
        active_power=0.5,
        reactive_power=1.1,
        R_th=1e-5,
        X_th=1e-5,
    )
    add_component!(sys, s2)
    g1 = ThermalStandard(
        name="init",
        available=true,
        status=false,
        bus=b2,
        active_power=0.1,
        reactive_power=0.1,
        rating=0.0,
        active_power_limits=(min=0.0, max=0.0),
        reactive_power_limits=nothing,
        ramp_limits=nothing,
        operation_cost=ThreePartCost(nothing),
        base_power=100.0,
        time_limits=nothing,
        prime_mover=PrimeMovers.OT,
        fuel=ThermalFuels.OTHER,
        services=Device[],
        dynamic_injector=nothing,
        ext=Dict{String, Any}(),
        time_series_container=InfrastructureSystems.TimeSeriesContainer(),
    )
    add_component!(sys, g1)

    @test run_powerflow!(sys, finite_diff=true)
    @test isapprox(
        get_active_power(get_component(Source, sys, "source_2")),
        0.5;
        atol=0.001,
    )
    @test isapprox(
        get_reactive_power(get_component(Source, sys, "source_2")),
        1.1;
        atol=0.001,
    )
end
