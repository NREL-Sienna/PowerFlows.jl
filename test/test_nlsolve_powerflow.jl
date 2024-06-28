@testset "NLsolve Power Flow 14-Bus testing" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    result_14 = [
        2.3255081760423684
        -0.15529254415401786
        0.4692141795968793
        -0.08704571860678871
        0.27136388547189727
        -0.22239752522709744
        1.014232467422514
        -0.17900878100582837
        1.0172382900002028
        -0.15297197376986682
        0.21603888720954267
        -0.2516366411330649
        1.0503438761049335
        -0.23128945395825606
        0.2453884476050094
        -0.23128945395825604
        1.0337109552890456
        -0.2588720100456853
        1.0325606430799592
        -0.26251860149671896
        1.0474837090281963
        -0.25914254888894855
        1.0535012072934
        -0.26648433793564014
        1.0471069158440354
        -0.2671769024662062
        1.0213119628726421
        -0.2803812119374241
    ]
    data = PowerFlows.PowerFlowData(ACPowerFlow(), sys; check_connectivity = true)
    #Compare results between finite diff methods and Jacobian method
    res1 = PowerFlows._solve_powerflow!(data, false)
    @test LinearAlgebra.norm(result_14 - res1.zero) <= 1e-6
    @test solve_ac_powerflow!(sys; method = :newton)

    # Test enforcing the reactive power Limits
    set_reactive_power!(get_component(PowerLoad, sys, "Bus4"), 0.0)
    data = PowerFlows.PowerFlowData(ACPowerFlow(), sys; check_connectivity = true)
    res2 = PowerFlows._solve_powerflow!(data, true)
    @test LinearAlgebra.norm(result_14 - res2.zero) >= 1e-6
    @test 1.08 <= res2.zero[15] <= 1.09
end

@testset "NLsolve Power Flow 14-Bus Line Configurations" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    base_res = solve_powerflow(ACPowerFlow(), sys)
    branch = first(PSY.get_components(Line, sys))
    dyn_branch = DynamicBranch(branch)
    add_component!(sys, dyn_branch)
    @test dyn_pf = solve_ac_powerflow!(sys)
    dyn_pf = solve_powerflow(ACPowerFlow(), sys)
    @test LinearAlgebra.norm(dyn_pf["bus_results"].Vm - base_res["bus_results"].Vm) <=
          1e-6

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    solve_ac_powerflow!(sys)
    @test PSY.get_active_power_flow(line) == 0.0
    test_bus = get_component(PSY.Bus, sys, "Bus 4")
    @test isapprox(PSY.get_magnitude(test_bus), 1.002; atol = 1e-3)

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    res = solve_powerflow(ACPowerFlow(), sys)
    @test res["flow_results"].P_from_to[4] == 0.0
    @test res["flow_results"].P_to_from[4] == 0.0
end

@testset "NLsolve Power Flow 3-Bus Fixed FixedAdmittance testing" begin
    p_gen_matpower_3bus = [20.3512373930753, 100.0, 100.0]
    q_gen_matpower_3bus = [45.516916781567232, 10.453799727283879, -31.992561631394636]
    sys_3bus = PSB.build_system(PSB.PSYTestSystems, "psse_3bus_gen_cls_sys")
    bus_103 = get_component(PSY.Bus, sys_3bus, "BUS 3")
    fix_shunt = PSY.FixedAdmittance("FixAdmBus3", true, bus_103, 0.0 + 0.2im)
    add_component!(sys_3bus, fix_shunt)
    df = solve_powerflow(ACPowerFlow(), sys_3bus)
    @test isapprox(df["bus_results"].P_gen, p_gen_matpower_3bus, atol = 1e-4)
    @test isapprox(df["bus_results"].Q_gen, q_gen_matpower_3bus, atol = 1e-4)
end

@testset "NLsolve Power Flow convergence fail testing" begin
    pf_sys5_re = PSB.build_system(PSB.PSITestSystems, "c_sys5_re"; add_forecasts = false)
    remove_component!(Line, pf_sys5_re, "1")
    remove_component!(Line, pf_sys5_re, "2")
    br = get_component(Line, pf_sys5_re, "6")
    PSY.set_x!(br, 20.0)
    PSY.set_r!(br, 2.0)

    # This is a negative test. The data passed for sys5_re is known to be infeasible.
    @test_logs(
        (:error, "The powerflow solver returned convergence = false"),
        match_mode = :any,
        @test !solve_ac_powerflow!(pf_sys5_re)
    )
end

@testset "Test 240 Case PSS/e results" begin
    file = joinpath(
        TEST_FILES_DIR,
        "test_data",
        "WECC240_v04_DPV_RE20_v33_6302_xfmr_DPbuscode_PFadjusted_V32_noRemoteVctrl.raw",
    )
    system = System(
        file;
        bus_name_formatter = x -> strip(string(x["name"])) * "-" * string(x["index"]),
        runchecks = false,
    )

    pf_bus_result_file = joinpath(TEST_FILES_DIR, "test_data", "pf_bus_results.csv")
    pf_gen_result_file = joinpath(TEST_FILES_DIR, "test_data", "pf_gen_results.csv")

    pf = solve_ac_powerflow!(system)
    @test pf
    pf_result_df = solve_powerflow(ACPowerFlow(), system)

    v_diff, angle_diff, number = psse_bus_results_compare(pf_bus_result_file, pf_result_df)
    p_diff, q_diff, names = psse_gen_results_compare(pf_gen_result_file, system)

    base_power = get_base_power(system)
    @test norm(v_diff, Inf) < DIFF_INF_TOLERANCE
    @test norm(v_diff, 2) / length(v_diff) < DIFF_L2_TOLERANCE
    @test norm(angle_diff, Inf) < DIFF_INF_TOLERANCE
    @test norm(angle_diff, 2) / length(angle_diff) < DIFF_L2_TOLERANCE
    @test norm(p_diff, Inf) < DIFF_INF_TOLERANCE * base_power
    @test norm(p_diff, 2) / length(p_diff) < DIFF_L2_TOLERANCE
    @test sum(q_diff) < DIFF_INF_TOLERANCE * base_power
    @test norm(q_diff, 2) / length(q_diff) < DIFF_L2_TOLERANCE
end

@testset "Multiple sources at ref" begin
    sys = System(100.0)
    b = ACBus(;
        number = 1,
        name = "01",
        bustype = ACBusTypes.REF,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b)

    #Test two sources with equal and opposite P and Q
    s1 = Source(;
        name = "source_1",
        available = true,
        bus = b,
        active_power = 0.5,
        reactive_power = 0.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s1)
    s2 = Source(;
        name = "source_2",
        available = true,
        bus = b,
        active_power = -0.5,
        reactive_power = -0.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s2)
    @test solve_ac_powerflow!(sys)

    #Create power mismatch, test for error
    set_active_power!(get_component(Source, sys, "source_1"), -0.4)
    @test_throws ErrorException(
        "Sources do not match P and/or Q requirements for reference bus.",
    ) solve_ac_powerflow!(sys)
end

@testset "AC PowerFlow with Multiple sources at PV" begin
    sys = System(100.0)
    b1 = ACBus(;
        number = 1,
        name = "01",
        bustype = ACBusTypes.REF,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b1)
    b2 = ACBus(;
        number = 2,
        name = "02",
        bustype = ACBusTypes.PV,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b2)
    a = Arc(; from = b1, to = b2)
    add_component!(sys, a)
    l = Line(;
        name = "l1",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = a,
        r = 1e-3,
        x = 1e-3,
        b = (from = 0.0, to = 0.0),
        rating = 0.0,
        angle_limits = (min = -pi / 2, max = pi / 2),
    )
    add_component!(sys, l)

    #Test two sources with equal and opposite P and Q
    s1 = Source(;
        name = "source_1",
        available = true,
        bus = b1,
        active_power = 0.5,
        reactive_power = 0.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s1)
    s2 = Source(;
        name = "source_2",
        available = true,
        bus = b2,
        active_power = 0.5,
        reactive_power = 1.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s2)
    s3 = Source(;
        name = "source_3",
        available = true,
        bus = b2,
        active_power = -0.5,
        reactive_power = -1.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s3)

    @test solve_ac_powerflow!(sys)

    #Create power mismatch, test for error
    set_reactive_power!(get_component(Source, sys, "source_3"), -0.5)
    @test_throws ErrorException("Sources do not match Q requirements for PV bus.") solve_ac_powerflow!(
        sys,
    )
end

@testset "AC PowerFlow Source + non-source at Ref" begin
    sys = System(100.0)
    b = ACBus(;
        number = 1,
        name = "01",
        bustype = ACBusTypes.REF,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b)

    #Test two sources with equal and opposite P and Q
    s1 = Source(;
        name = "source_1",
        available = true,
        bus = b,
        active_power = 0.5,
        reactive_power = 0.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s1)
    g1 = ThermalStandard(;
        name = "init",
        available = true,
        status = false,
        bus = b,
        active_power = 0.1,
        reactive_power = 0.1,
        rating = 0.0,
        active_power_limits = (min = 0.0, max = 0.0),
        reactive_power_limits = nothing,
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
    add_component!(sys, g1)

    @test solve_ac_powerflow!(sys)
    @test isapprox(
        get_active_power(get_component(Source, sys, "source_1")),
        0.5;
        atol = 0.001,
    )
    @test isapprox(
        get_reactive_power(get_component(Source, sys, "source_1")),
        0.1;
        atol = 0.001,
    )
end

@testset "AC PowerFlow Source + non-source at PV" begin
    sys = System(100.0)
    b1 = ACBus(;
        number = 1,
        name = "01",
        bustype = ACBusTypes.REF,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b1)
    b2 = ACBus(;
        number = 2,
        name = "02",
        bustype = ACBusTypes.PV,
        angle = 0.0,
        magnitude = 1.1,
        voltage_limits = (0.0, 2.0),
        base_voltage = 230,
    )
    add_component!(sys, b2)
    a = Arc(; from = b1, to = b2)
    add_component!(sys, a)
    l = Line(;
        name = "l1",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = a,
        r = 1e-3,
        x = 1e-3,
        b = (from = 0.0, to = 0.0),
        rating = 0.0,
        angle_limits = (min = -pi / 2, max = pi / 2),
    )
    add_component!(sys, l)

    #Test two sources with equal and opposite P and Q
    s1 = Source(;
        name = "source_1",
        available = true,
        bus = b1,
        active_power = 0.5,
        reactive_power = 0.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s1)
    s2 = Source(;
        name = "source_2",
        available = true,
        bus = b2,
        active_power = 0.5,
        reactive_power = 1.1,
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, s2)
    g1 = ThermalStandard(;
        name = "init",
        available = true,
        status = false,
        bus = b2,
        active_power = 0.1,
        reactive_power = 0.1,
        rating = 0.0,
        active_power_limits = (min = 0.0, max = 0.0),
        reactive_power_limits = nothing,
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
    add_component!(sys, g1)

    @test solve_ac_powerflow!(sys)
    @test isapprox(
        get_active_power(get_component(Source, sys, "source_2")),
        0.5;
        atol = 0.001,
    )
    @test isapprox(
        get_reactive_power(get_component(Source, sys, "source_2")),
        1.1;
        atol = 0.001,
    )
end
