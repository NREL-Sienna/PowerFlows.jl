@testset "AC Power Flow 14-Bus testing" for ACSolver in AC_SOLVERS_TO_TEST
    result_14 = [
        2.3255081760423684
        -0.15529254415401786
        0.4692141795968793 - 0.127  # Q_gen - Q_load
        -0.08704571860678871
        0.27136388547189727 - 0.19  # Q_gen - Q_load
        -0.22239752522709744
        1.014232467422514
        -0.17900878100582837
        1.0172382900002028
        -0.15297197376986682
        0.21603888720954267 - 0.075  # Q_gen - Q_load
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

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    set_units_base_system!(sys, UnitSystem.SYSTEM_BASE)
    pf = ACPowerFlow{ACSolver}()
    data = PowerFlows.PowerFlowData(pf, sys; check_connectivity = true)
    #Compare results between finite diff methods and Jacobian method
    converged1 = PowerFlows._ac_powerflow(data, pf, 1)
    x1 = _calc_x(data, 1)
    @test LinearAlgebra.norm(result_14 - x1, Inf) <= 1e-6

    # Test that solve_powerflow! succeeds
    solved1 = deepcopy(sys)
    @test solve_powerflow!(pf, solved1)

    # Test that passing check_reactive_power_limits=false is the default and violates limits
    solved2 = deepcopy(sys)
    @test solve_powerflow!(pf, solved2; check_reactive_power_limits = false)
    @test IS.compare_values(solved1, solved2)
    @test get_reactive_power(get_component(ThermalStandard, solved2, "Bus8")) >
          get_reactive_power_limits(get_component(ThermalStandard, solved2, "Bus8")).max

    # Test that passing check_reactive_power_limits=true fixes that
    solved3 = deepcopy(sys)
    @test solve_powerflow!(pf, solved3; check_reactive_power_limits = true)
    @test get_reactive_power(get_component(ThermalStandard, solved3, "Bus8")) <=
          get_reactive_power_limits(get_component(ThermalStandard, solved3, "Bus8")).max

    # Test Newton method
    @test solve_powerflow!(pf, deepcopy(sys))

    # Test enforcing the reactive power limits in closer detail
    set_reactive_power!(get_component(PowerLoad, sys, "Bus4"), 0.0)
    data = PowerFlows.PowerFlowData(pf, sys; check_connectivity = true)
    converged2 = PowerFlows._ac_powerflow(data, pf, 1; check_reactive_power_limits = true)
    x2 = _calc_x(data, 1)
    @test LinearAlgebra.norm(result_14 - x2, Inf) >= 1e-6
    @test 1.08 <= x2[15] <= 1.09
end

@testset "AC Power Flow 14-Bus Line Configurations" for ACSolver in AC_SOLVERS_TO_TEST
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    pf = ACPowerFlow{ACSolver}()
    base_res = solve_powerflow(pf, sys)
    branch = first(PSY.get_components(Line, sys))
    dyn_branch = DynamicBranch(branch)
    add_component!(sys, dyn_branch)
    @test dyn_pf = solve_powerflow!(pf, sys)
    dyn_pf = solve_powerflow(pf, sys)
    @test LinearAlgebra.norm(dyn_pf["bus_results"].Vm - base_res["bus_results"].Vm, Inf) <=
          1e-6

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    solve_powerflow!(pf, sys)
    @test PSY.get_active_power_flow(line) == 0.0
    test_bus = get_component(PSY.Bus, sys, "Bus 4")
    @test isapprox(PSY.get_magnitude(test_bus), 1.002; atol = 1e-3, rtol = 0)

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    res = solve_powerflow(pf, sys)
    @test res["flow_results"].P_from_to[4] == 0.0
    @test res["flow_results"].P_to_from[4] == 0.0
end

@testset "AC Power Flow 3-Bus Fixed FixedAdmittance testing" for ACSolver in
                                                                 AC_SOLVERS_TO_TEST
    p_gen_matpower_3bus = [20.3512373930753, 100.0, 100.0]
    q_gen_matpower_3bus = [45.516916781567232, 10.453799727283879, -31.992561631394636]
    sys_3bus = PSB.build_system(PSB.PSYTestSystems, "psse_3bus_gen_cls_sys")
    bus_103 = get_component(PSY.Bus, sys_3bus, "BUS 3")
    fix_shunt = PSY.FixedAdmittance("FixAdmBus3", true, bus_103, 0.0 + 0.2im)
    add_component!(sys_3bus, fix_shunt)
    pf = ACPowerFlow{ACSolver}()
    df = solve_powerflow(pf, sys_3bus)
    @test isapprox(df["bus_results"].P_gen, p_gen_matpower_3bus, atol = 1e-4)
    @test isapprox(df["bus_results"].Q_gen, q_gen_matpower_3bus, atol = 1e-4)
end

@testset "AC Power Flow convergence fail testing" for ACSolver in AC_SOLVERS_TO_TEST
    pf_sys5_re = PSB.build_system(PSB.PSITestSystems, "c_sys5_re"; add_forecasts = false)
    remove_component!(Line, pf_sys5_re, "1")
    remove_component!(Line, pf_sys5_re, "2")
    br = get_component(Line, pf_sys5_re, "6")
    PSY.set_x!(br, 20.0)
    PSY.set_r!(br, 2.0)

    pf = ACPowerFlow{ACSolver}()

    # This is a negative test. The data passed for sys5_re is known to be infeasible.
    @test_logs(
        (:error, "The powerflow solver returned convergence = false"),
        match_mode = :any,
        @test !solve_powerflow!(pf, pf_sys5_re)
    )
end

@testset "AC Test 240 Case PSS/e results" for ACSolver in AC_SOLVERS_TO_TEST
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

    pf = ACPowerFlow{ACSolver}()

    pf1 = solve_powerflow!(pf, system)
    @test pf1
    pf_result_df = solve_powerflow(pf, system)

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

@testset "AC Multiple sources at ref" for ACSolver in AC_SOLVERS_TO_TEST
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
    pf = ACPowerFlow{ACSolver}()
    @test solve_powerflow!(pf, sys)

    #Create power mismatch, test for error
    set_active_power!(get_component(Source, sys, "source_1"), -0.4)
    @test_throws ErrorException(
        "Sources do not match P and/or Q requirements for reference bus.",
    ) solve_powerflow!(pf, sys)
end

@testset "AC PowerFlow with Multiple sources at PV" for ACSolver in AC_SOLVERS_TO_TEST
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

    pf = ACPowerFlow{ACSolver}()

    @test solve_powerflow!(pf, sys)

    #Create power mismatch, test for error
    set_reactive_power!(get_component(Source, sys, "source_3"), -0.5)
    @test_throws ErrorException("Sources do not match Q requirements for PV bus.") solve_powerflow!(
        pf,
        sys,
    )
end

@testset "AC PowerFlow Source + non-source at Ref" for ACSolver in AC_SOLVERS_TO_TEST
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

    pf = ACPowerFlow{ACSolver}()

    @test solve_powerflow!(pf, sys)
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

@testset "AC PowerFlow Source + non-source at PV" for ACSolver in AC_SOLVERS_TO_TEST
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

    pf = ACPowerFlow{ACSolver}()

    @test solve_powerflow!(pf, sys)
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

# in this test, the following aspects are checked:
# 1. The results of the power flow are consistent for the KLU and Hybrid solvers
# 2. The results of the power flow are consistent for the KLU solver and the legacy implementation
# 3. The Jacobian matrix is the same for the KLU solver and the legacy implementation
@testset "Compare larger grid results KLU vs Hybrid" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")

    pf_default = ACPowerFlow()
    pf_lu = ACPowerFlow(LUACPowerFlow)
    pf_newton = ACPowerFlow(NewtonRaphsonACPowerFlow)

    PSY.set_units_base_system!(sys, "SYSTEM_BASE")
    data = PowerFlowData(
        pf_default,
        sys;
        check_connectivity = true)

    time_step = 1

    res_default = solve_powerflow(pf_default, sys)  # must be the same as KLU
    res_lu = solve_powerflow(pf_lu, sys)
    res_newton = solve_powerflow(pf_newton, sys)

    @test all(
        isapprox.(
            res_lu["bus_results"][!, :Vm],
            res_default["bus_results"][!, :Vm],
            rtol = 0,
            atol = 1e-12,
        ),
    )
    @test all(
        isapprox.(
            res_lu["bus_results"][!, :θ],
            res_default["bus_results"][!, :θ],
            rtol = 0,
            atol = 1e-12,
        ),
    )

    @test all(
        isapprox.(
            res_lu["bus_results"][!, :Vm],
            res_newton["bus_results"][!, :Vm],
            rtol = 0,
            atol = 1e-12,
        ),
    )
    @test all(
        isapprox.(
            res_lu["bus_results"][!, :θ],
            res_newton["bus_results"][!, :θ],
            rtol = 0,
            atol = 1e-12,
        ),
    )
end

@testset "Test loss factors for larger grid" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")

    pf_lu = ACPowerFlow(LUACPowerFlow)
    pf_lu_lf = ACPowerFlow(LUACPowerFlow; calculate_loss_factors = true)
    pf_newton = ACPowerFlow(NewtonRaphsonACPowerFlow; calculate_loss_factors = true)

    data_lu = PowerFlowData(
        pf_lu_lf,
        sys;
        check_connectivity = true)

    data_newton = PowerFlowData(
        pf_newton,
        sys;
        check_connectivity = true)

    data_brute_force = PowerFlowData(
        pf_newton,
        sys;
        check_connectivity = true)

    time_step = 1

    solve_powerflow!(data_lu; pf = pf_lu)
    solve_powerflow!(data_newton; pf = pf_newton)

    @test all(
        isapprox.(
            data_lu.loss_factors,
            data_newton.loss_factors,
            rtol = 0,
            atol = 1e-9,
        ),
    )

    bf_loss_factors =
        penalty_factors_brute_force(data_brute_force, pf_newton)
    @test all(isapprox.(
        data_newton.loss_factors,
        bf_loss_factors,
        rtol = 0,
        atol = 1e-4,
    ))
end

@testset "AC PF with distributed slack" for (grid_lib, grid_name) in [
    (PSB.PSITestSystems, "c_sys14"),
    (PSB.MatpowerTestSystems, "matpower_case30_sys"),
]
    function _get_spf_dict(slack_participation_factors)
        generator_slack_participation_factors = Dict{String, Float64}()
        for (b, spf) in enumerate(slack_participation_factors)
            get_bustype(
                get_bus(sys, bus_numbers[b])
            ) ==
            ACBusTypes.PQ && continue
            gens = get_components(
                x -> get_number(get_bus(x)) == bus_numbers[b],
                ThermalStandard,
                sys,
            )
            isempty(gens) && continue
            gens = collect(gens)
            for g in gens
                generator_slack_participation_factors[get_name(g)] = spf /
                                                                     length(gens)
            end
        end
        return generator_slack_participation_factors
    end

    sys = PSB.build_system(grid_lib, grid_name)

    # add a duplicate generator to a PV bus to make sure the method works for such set-ups
    g1 = first(
        get_components(x -> get_bustype(get_bus(x)) == ACBusTypes.PV, ThermalStandard, sys),
    )

    g2 = ThermalStandard(;
        name = "Duplicate",
        available = true,
        status = true,
        bus = get_bus(g1),
        active_power = 0.1,
        reactive_power = 0.1,
        rating = 1.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -1.0, max = 1.0),
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
    add_component!(sys, g2)

    bus_numbers = get_bus_numbers(sys)

    ref_n = []
    pv_n = []
    for (i, bn) in enumerate(bus_numbers)
        isempty(get_components(x -> get_number(get_bus(x)) == bn, ThermalStandard, sys)) &&
            continue
        b = only(get_components(x -> get_number(x) == bn, ACBus, sys))
        get_bustype(b) == ACBusTypes.REF && (push!(ref_n, i))
        get_bustype(b) == ACBusTypes.PV && (push!(pv_n, i))
    end

    # make sure we have active power imbalance in the starting grid
    g = first(
        get_components(
            x -> get_bustype(get_bus(x)) == ACBusTypes.REF,
            ThermalStandard,
            sys,
        ),
    )
    with_units_base(sys, UnitSystem.NATURAL_UNITS) do
        set_active_power!(g, 20.0)
    end

    pf = ACPowerFlow()
    data = PowerFlowData(pf, sys)
    original_bus_power, original_gen_power = _system_generation_power(sys, bus_numbers)
    data_original_bus_power = copy(data.bus_activepower_injection[:, 1])
    res1 = solve_powerflow(pf, sys)

    slack_participation_factors = zeros(Float64, length(bus_numbers))
    slack_participation_factors[ref_n] .= 1.0

    pf2 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(slack_participation_factors),
    )
    res2 = solve_powerflow(pf2, sys)

    # basic test: if we pass the same slack participation factors as the default ones, the results
    # should be the same
    @test isapprox(res1["bus_results"].Vm, res2["bus_results"].Vm, atol = 1e-6, rtol = 0)
    @test isapprox(res1["bus_results"].θ, res2["bus_results"].θ, atol = 1e-6, rtol = 0)

    @test _check_ds_pf(
        pf2,
        sys,
        slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # now test with REF and one PV bus having slack participation factors of 1.0
    slack_participation_factors[pv_n[1]] = 1.0
    pf3 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(slack_participation_factors))

    @test _check_ds_pf(
        pf3,
        sys,
        slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # now test with all REF and PV buses having equal slack participation factors of 1.0
    slack_participation_factors[pv_n] .= 1.0
    pf4 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(slack_participation_factors))

    @test _check_ds_pf(
        pf4,
        sys,
        slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # Now set the slack participation factor to 0.0 for the REF bus
    slack_participation_factors[ref_n] .= 0.0
    pf5 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(slack_participation_factors))

    @test _check_ds_pf(
        pf5,
        sys,
        slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # now check the formula of the distribution of slack provision for different factors
    slack_participation_factors[ref_n] .= 2.5
    slack_participation_factors[pv_n] .= pv_n
    pf6 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(slack_participation_factors))

    @test _check_ds_pf(
        pf6,
        sys,
        slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # TODO: add a test when two generators at the same bus have unequal slack participation factors
    # TODO: add test for multiperiod with just one time step and one Dict input in the vector
    # TODO: add test for multiperiod with just one Dict input of slack factors for different time steps
    # TODO: add test for multiperiod with different Dict inputs for different time steps
end
