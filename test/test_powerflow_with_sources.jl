@testset "Powerflow with multiple sources at bus" begin
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
