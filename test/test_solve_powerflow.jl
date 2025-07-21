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
    data = PowerFlows.PowerFlowData(
        pf,
        sys;
        check_connectivity = true,
        correct_bustypes = true,
    )
    #Compare results between finite diff methods and Jacobian method
    converged1 = PowerFlows._ac_powerflow(data, pf, 1)
    x1 = _calc_x(data, 1)
    @test LinearAlgebra.norm(result_14 - x1, Inf) <= 1e-6 # <- this fails likely due to the change of the B allocation
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
    data = PowerFlows.PowerFlowData(
        pf,
        sys;
        check_connectivity = true,
        correct_bustypes = true,
    )
    converged2 = PowerFlows._ac_powerflow(data, pf, 1; check_reactive_power_limits = true)
    x2 = _calc_x(data, 1)
    @test LinearAlgebra.norm(result_14 - x2, Inf) >= 1e-6
    @test 1.08 <= x2[15] <= 1.09
end

@testset "AC Power Flow 14-Bus Line Configurations" for ACSolver in AC_SOLVERS_TO_TEST
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    pf = ACPowerFlow{ACSolver}()
    base_res = solve_powerflow(pf, sys; correct_bustypes = true)
    branch = first(PSY.get_components(Line, sys))
    dyn_branch = DynamicBranch(branch)
    add_component!(sys, dyn_branch)
    @test dyn_pf = solve_powerflow!(pf, sys; correct_bustypes = true)
    dyn_pf = solve_powerflow(pf, sys; correct_bustypes = true)
    @test LinearAlgebra.norm(dyn_pf["bus_results"].Vm - base_res["bus_results"].Vm, Inf) <=
          1e-6

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    solve_powerflow!(pf, sys; correct_bustypes = true)
    @test PSY.get_active_power_flow(line) == 0.0
    test_bus = get_component(PSY.ACBus, sys, "Bus 4")
    @test isapprox(PSY.get_magnitude(test_bus), 1.002; atol = 1e-3, rtol = 0)

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    res = solve_powerflow(pf, sys; correct_bustypes = true)
    @test res["flow_results"].P_from_to[4] == 0.0
    @test res["flow_results"].P_to_from[4] == 0.0
end

@testset "AC Power Flow 3-Bus Fixed FixedAdmittance testing" for ACSolver in
                                                                 AC_SOLVERS_TO_TEST
    p_gen_matpower_3bus = [20.3512373930753, 100.0, 100.0]
    q_gen_matpower_3bus = [45.516916781567232, 10.453799727283879, -31.992561631394636]
    sys_3bus = PSB.build_system(PSB.PSYTestSystems, "psse_3bus_gen_cls_sys")
    bus_103 = get_component(PSY.ACBus, sys_3bus, "BUS 3")
    fix_shunt = PSY.FixedAdmittance("FixAdmBus3", true, bus_103, 0.0 + 0.2im)
    add_component!(sys_3bus, fix_shunt)
    pf = ACPowerFlow{ACSolver}()
    df = solve_powerflow(pf, sys_3bus; correct_bustypes = true)
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

# FIXME currently errors: write_powerflow_solution! relies on all PV buses in
# the system having available generators.
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

    pf1 = solve_powerflow!(pf, system; correct_bustypes = true)
    @test pf1
    pf_result_df = solve_powerflow(pf, system; correct_bustypes = true)

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
    b = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)

    #Test two sources with equal and opposite P and Q
    s1 = _add_simple_source!(sys, b, 0.5, 0.1)
    s2 = _add_simple_source!(sys, b, -0.5, -0.1)

    pf = ACPowerFlow{ACSolver}()
    @test solve_powerflow!(pf, sys; correct_bustypes = true)

    #Create power mismatch, test for error
    set_active_power!(s1, -0.4)
    @test_throws ErrorException(
        "Sources do not match P and/or Q requirements for reference bus.",
    ) solve_powerflow!(pf, sys)
end

@testset "AC PowerFlow with Multiple sources at PV" for ACSolver in AC_SOLVERS_TO_TEST
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 2, ACBusTypes.PV, 230, 1.1, 0.0)

    l = _add_simple_line!(sys, b1, b2, 1e-3, 1e-3, 0.0)

    #Test two sources with equal and opposite P and Q
    s1 = _add_simple_source!(sys, b1, 0.5, 0.1)
    s2 = _add_simple_source!(sys, b2, 0.5, 1.1)
    s3 = _add_simple_source!(sys, b2, -0.5, -1.1)

    pf = ACPowerFlow{ACSolver}()

    @test solve_powerflow!(pf, sys; correct_bustypes = true)

    #Create power mismatch, test for error
    set_reactive_power!(s3, -0.5)
    @test_throws ErrorException("Sources do not match Q requirements for PV bus.") solve_powerflow!(
        pf,
        sys,
        correct_bustypes = true,
    )
end

@testset "AC PowerFlow Source + non-source at Ref" for ACSolver in AC_SOLVERS_TO_TEST
    sys = System(100.0)

    b = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)

    #Test two sources with equal and opposite P and Q
    s1 = _add_simple_source!(sys, b, 0.5, 0.1)
    g1 = _add_simple_thermal_standard!(sys, b, 0.1, 0.1)

    pf = ACPowerFlow{ACSolver}()

    @test solve_powerflow!(pf, sys)
    @test isapprox(get_active_power(s1), 0.5; atol = 0.001)
    @test isapprox(get_reactive_power(s1), 0.1; atol = 0.001)
end

@testset "AC PowerFlow Source + non-source at PV" for ACSolver in AC_SOLVERS_TO_TEST
    sys = System(100.0)

    b1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 2, ACBusTypes.PV, 230, 1.1, 0.0)
    l = _add_simple_line!(sys, b1, b2, 1e-3, 1e-3, 0.0)

    #Test two sources with equal and opposite P and Q
    s1 = _add_simple_source!(sys, b1, 0.5, 0.1)
    s2 = _add_simple_source!(sys, b2, 0.5, 1.1)
    g1 = _add_simple_thermal_standard!(sys, b2, 0.1, 0.1)

    pf = ACPowerFlow{ACSolver}()

    @test solve_powerflow!(pf, sys; correct_bustypes = true)
    @test isapprox(get_active_power(s2), 0.5; atol = 0.001)
    @test isapprox(get_reactive_power(s2), 1.1; atol = 0.001)
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
        check_connectivity = true,
        correct_bustypes = true)

    time_step = 1

    res_default = solve_powerflow(pf_default, sys; correct_bustypes = true)  # must be the same as KLU
    res_lu = solve_powerflow(pf_lu, sys; correct_bustypes = true)
    res_newton = solve_powerflow(pf_newton, sys; correct_bustypes = true)

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
    pf_lu_lf = ACPowerFlow(
        LUACPowerFlow;
        calculate_loss_factors = true,
        calculate_voltage_stability_factors = true,
    )
    pf_newton = ACPowerFlow(
        NewtonRaphsonACPowerFlow;
        calculate_loss_factors = true,
        calculate_voltage_stability_factors = true,
    )

    data_lu = PowerFlowData(
        pf_lu_lf,
        sys;
        check_connectivity = true,
        correct_bustypes = true)

    data_newton = PowerFlowData(
        pf_newton,
        sys;
        check_connectivity = true,
        correct_bustypes = true)

    data_brute_force = PowerFlowData(
        pf_newton,
        sys;
        check_connectivity = true,
        correct_bustypes = true)

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

    @test all(
        isapprox.(
            data_lu.voltage_stability_factors,
            data_newton.voltage_stability_factors,
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

@testset "voltage_stability_factors" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    pf_lu = ACPowerFlow(LUACPowerFlow; calculate_voltage_stability_factors = true)
    pf_newton =
        ACPowerFlow(NewtonRaphsonACPowerFlow; calculate_voltage_stability_factors = true)
    data_lu = PowerFlowData(
        pf_lu,
        sys;
        check_connectivity = true,
        correct_bustypes = true,
    )
    data_newton = PowerFlowData(
        pf_newton,
        sys;
        check_connectivity = true,
        correct_bustypes = true,
    )
    time_step = 1
    solve_powerflow!(data_lu; pf = pf_lu)
    solve_powerflow!(data_newton; pf = pf_newton)
    @test all(
        isapprox.(
            data_lu.voltage_stability_factors,
            data_newton.voltage_stability_factors,
            rtol = 0,
            atol = 1e-6,
        ),
    )
    ref, pv, pq = PowerFlows.bus_type_idx(data_lu, time_step)
    pvpq = [pv; pq]
    npvpq = length(pvpq)
    V = data_lu.bus_magnitude[:, time_step] .* exp.(1im * data_lu.bus_angles[:, time_step])
    dSbus_dVa, dSbus_dVm = _legacy_dSbus_dV(V, data_lu.power_network_matrix.data)
    J = _legacy_J(dSbus_dVa, dSbus_dVm, pvpq, pq)
    Gs =
        J[(npvpq + 1):end, (npvpq + 1):end] -
        J[(npvpq + 1):end, 1:npvpq] * inv(collect(J[1:npvpq, 1:npvpq])) *
        J[1:npvpq, (npvpq + 1):end]
    u_1, (σ_1,), v_1, _ = PROPACK.tsvd_irl(Gs; smallest = true, k = 1)
    σ, u, v = PowerFlows._singular_value_decomposition(J, npvpq)

    @assert isapprox(σ_1, σ, atol = 1e-6)
    # the sign does not matter
    @assert isapprox(sign(first(u_1)) * u_1, u, atol = 1e-4)
    @assert isapprox(sign(first(v_1)) * v_1, v, atol = 1e-4)
end

@testset "AC PF with distributed slack" for (grid_lib, grid_name) in [
        (PSB.PSITestSystems, "c_sys14"),
        (PSB.MatpowerTestSystems, "matpower_case30_sys"),
    ], ACSolver in (NewtonRaphsonACPowerFlow, TrustRegionACPowerFlow)
    function _get_spf_dict(bus_slack_participation_factors)
        generator_slack_participation_factors = Dict{Tuple{DataType, String}, Float64}()
        for (b, spf) in enumerate(bus_slack_participation_factors)
            get_bustype(get_bus(sys, bus_numbers[b])) == ACBusTypes.PQ && continue
            gens = get_components(
                x -> get_number(get_bus(x)) == bus_numbers[b],
                ThermalStandard,
                sys,
            )
            isempty(gens) && continue
            gens = collect(gens)
            for g in gens
                generator_slack_participation_factors[(ThermalStandard, get_name(g))] =
                    spf / length(gens)
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
        bus_type = get_bustype(b)
        bus_type == ACBusTypes.REF && (push!(ref_n, i))
        bus_type == ACBusTypes.PV && (push!(pv_n, i))
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
    data = PowerFlowData(pf, sys; correct_bustypes = true)
    original_bus_power, original_gen_power = _system_generation_power(sys, bus_numbers)
    data_original_bus_power = copy(data.bus_activepower_injection[:, 1])
    res1 = solve_powerflow(pf, sys; correct_bustypes = true)

    bus_slack_participation_factors = zeros(Float64, length(bus_numbers))
    bus_slack_participation_factors[ref_n] .= 1.0

    pf2 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(
            bus_slack_participation_factors,
        ),
    )
    res2 = solve_powerflow(pf2, sys; correct_bustypes = true)

    # basic test: if we pass the same slack participation factors as the default ones, the results
    # should be the same
    @test isapprox(res1["bus_results"].Vm, res2["bus_results"].Vm, atol = 1e-6, rtol = 0)
    @test isapprox(res1["bus_results"].θ, res2["bus_results"].θ, atol = 1e-6, rtol = 0)

    _check_ds_pf(
        pf2,
        sys,
        bus_slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # now test with REF and one PV bus having slack participation factors of 1.0
    bus_slack_participation_factors[pv_n[1]] = 1.0
    pf3 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(
            bus_slack_participation_factors,
        ),
    )

    _check_ds_pf(
        pf3,
        sys,
        bus_slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # now test with all REF and PV buses having equal slack participation factors of 1.0
    bus_slack_participation_factors[pv_n] .= 1.0
    pf4 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(
            bus_slack_participation_factors,
        ),
    )

    _check_ds_pf(
        pf4,
        sys,
        bus_slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # Now set the slack participation factor to 0.0 for the REF bus
    bus_slack_participation_factors[ref_n] .= 0.0
    pf5 = ACPowerFlow(;
        generator_slack_participation_factors = _get_spf_dict(
            bus_slack_participation_factors,
        ),
    )

    _check_ds_pf(
        pf5,
        sys,
        bus_slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )

    # now check the formula of the distribution of slack provision for different factors
    bus_slack_participation_factors[ref_n] .= 2.5
    bus_slack_participation_factors[pv_n] .= pv_n
    pf6 = ACPowerFlow(;
        generator_slack_participation_factors = [
            _get_spf_dict(
                bus_slack_participation_factors,
            ),
        ],
    )  # [] to test this input variant

    _check_ds_pf(
        pf6,
        sys,
        bus_slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power,
    )
end

@testset "AC PF DS power redistribution" for ACSolver in (
    NewtonRaphsonACPowerFlow,
    TrustRegionACPowerFlow,
)
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 2, ACBusTypes.PV, 230, 1.1, 0.0)
    l = _add_simple_line!(sys, b1, b2, 1e-3, 1e-3, 0.0)

    ps = -0.5
    s1 = _add_simple_source!(sys, b1, ps, 0.1)

    p1 = 0.1
    g1 = _add_simple_thermal_standard!(sys, b2, p1, 0.1)

    p2 = 0.2
    g2 = _add_simple_thermal_standard!(sys, b2, p2, 0.1)

    reset_p() =
        for (c, p) in zip((s1, g1, g2), (ps, p1, p2))
            set_active_power!(c, p)
        end

    gspf = Dict(
        (Source, get_name(s1)) => 1.0,
        (ThermalStandard, get_name(g1)) => 0.0,
        (ThermalStandard, get_name(g2)) => 0.0,
    )
    pf = ACPowerFlow(; generator_slack_participation_factors = gspf)
    solve_powerflow!(pf, sys; correct_bustypes = true)
    @test isapprox(get_active_power(g1), p1; atol = 1e-6, rtol = 0)
    @test isapprox(get_active_power(g2), p2; atol = 1e-6, rtol = 0)
    reset_p()

    gspf = Dict(
        (Source, get_name(s1)) => 0.0,
        (ThermalStandard, get_name(g1)) => 0.5,
        (ThermalStandard, get_name(g2)) => 0.5,
    )
    pf = ACPowerFlow(; generator_slack_participation_factors = gspf)
    solve_powerflow!(pf, sys; correct_bustypes = true)
    @test isapprox(get_active_power(s1), ps; atol = 1e-6, rtol = 0)
    @test isapprox(
        get_active_power(g1) - p1,
        get_active_power(g2) - p2;
        atol = 1e-6,
        rtol = 0,
    )
    reset_p()

    gspf =
        Dict((ThermalStandard, get_name(g1)) => 0.0, (ThermalStandard, get_name(g2)) => 1.0)
    pf = ACPowerFlow(; generator_slack_participation_factors = gspf)
    solve_powerflow!(pf, sys)
    @test isapprox(get_active_power(s1), ps; atol = 1e-6, rtol = 0)
    @test isapprox(get_active_power(g1), p1; atol = 1e-6, rtol = 0)
    @test isapprox(
        -get_active_power(g2),
        get_active_power(g1) + get_active_power(s1);
        atol = 1e-3,  # losses don't allow lower tolerance
        rtol = 0,
    )
    reset_p()

    gspf = Dict(
        (Source, get_name(s1)) => 0.1,
        (ThermalStandard, get_name(g1)) => 0.2,
        (ThermalStandard, get_name(g2)) => 0.4,
    )
    pf = ACPowerFlow(; generator_slack_participation_factors = gspf)
    solve_powerflow!(pf, sys)
    total_slack_power =
        -(get_active_power(s1) - ps) + get_active_power(g1) - p1 + get_active_power(g2) - p2
    @test isapprox(
        get_active_power(s1) - ps,
        total_slack_power * 0.2;
        atol = 1e-6,
        rtol = 0,
    )
    @test isapprox(
        get_active_power(g1) - p1,
        total_slack_power * 0.4;
        atol = 1e-6,
        rtol = 0,
    )
    @test isapprox(
        get_active_power(g2) - p2,
        total_slack_power * 0.8;
        atol = 1e-3,
        rtol = 0,
    )
    reset_p()
end

@testset "AC PF DS with two connected components" for mode in (:same, :random),
    gen_mode in (:gen, :source),
    ACSolver in (NewtonRaphsonACPowerFlow, TrustRegionACPowerFlow)
    # here we build two identical grids in one system
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 8, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 4, ACBusTypes.PV, 230, 1.1, 0.0)
    b3 = _add_simple_bus!(sys, 3, ACBusTypes.PQ, 230, 1.1, 0.0)

    b4 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b5 = _add_simple_bus!(sys, 5, ACBusTypes.PV, 230, 1.1, 0.0)
    b6 = _add_simple_bus!(sys, 2, ACBusTypes.PQ, 230, 1.1, 0.0)

    s1 = _add_simple_source!(sys, b1, 0.0, 0.0)
    s2 = _add_simple_source!(sys, b4, 0.0, 0.0)

    g1 = if gen_mode == :gen
        _add_simple_thermal_standard!(sys, b2, 0.0, 0.0)
    else
        _add_simple_source!(sys, b2, 0.0, 0.0)
    end
    g2 = if gen_mode == :gen
        _add_simple_thermal_standard!(sys, b5, 0.0, 0.0)
    else
        _add_simple_source!(sys, b5, 0.0, 0.0)
    end

    ld1 = _add_simple_load!(sys, b3, 6, 2)
    ld2 = _add_simple_load!(sys, b6, 6, 2)

    l1 = _add_simple_line!(sys, b1, b2, 1e-3, 1e-3, 0.0)
    l2 = _add_simple_line!(sys, b2, b3, 1e-3, 1e-3, 0.0)

    l3 = _add_simple_line!(sys, b4, b5, 1e-3, 1e-3, 0.0)
    l4 = _add_simple_line!(sys, b5, b6, 1e-3, 1e-3, 0.0)

    devices = collect(get_components(StaticInjection, sys))

    if mode == :same
        factors = ones(Float64, length(devices))
    elseif mode == :random
        Random.seed!(0)
        factors = abs.(randn(Float64, length(devices)))
    else
        error("Unknown mode: $mode")
    end

    generator_slack_participation_factors =
        Dict((typeof(x), get_name(x)) => f for (x, f) in zip(devices, factors))

    pf = ACPowerFlow(ACSolver;
        generator_slack_participation_factors = generator_slack_participation_factors,
    )

    data = PowerFlowData(pf, sys; check_connectivity = false, correct_bustypes = true)
    data_original_bus_power = copy(data.bus_activepower_injection[:, 1])
    bus_numbers = get_bus_numbers(sys)
    original_bus_power, original_gen_power = _system_generation_power(sys, bus_numbers)

    solve_powerflow!(pf, sys; check_connectivity = false, correct_bustypes = true)

    if mode == :same
        @test isapprox(get_active_power(s1), get_active_power(g1), rtol = 0, atol = 1e-6)
        @test isapprox(get_active_power(s2), get_active_power(g2), rtol = 0, atol = 1e-6)

        @test isapprox(get_active_power(s1), get_active_power(s2), rtol = 0, atol = 1e-6)
        @test isapprox(get_active_power(g1), get_active_power(g2), rtol = 0, atol = 1e-6)
    end

    # tol of 1e-3 due to losses
    @test isapprox(
        get_active_power(ld1),
        get_active_power(s1) + get_active_power(g1),
        rtol = 0,
        atol = 1e-3,
    )
    @test isapprox(
        get_active_power(ld2),
        get_active_power(s2) + get_active_power(g2),
        rtol = 0,
        atol = 1e-3,
    )

    solve_powerflow!(data; pf = pf)

    @test isapprox(
        data.bus_activepower_injection[data.bus_lookup[get_number(b1)], 1],
        get_active_power(s1),
        rtol = 0,
        atol = 1e-6,
    )
    @test isapprox(
        data.bus_activepower_injection[data.bus_lookup[get_number(b2)], 1],
        get_active_power(g1),
        rtol = 0,
        atol = 1e-6,
    )
    @test isapprox(
        data.bus_activepower_injection[data.bus_lookup[get_number(b4)], 1],
        get_active_power(s2),
        rtol = 0,
        atol = 1e-6,
    )
    @test isapprox(
        data.bus_activepower_injection[data.bus_lookup[get_number(b5)], 1],
        get_active_power(g2),
        rtol = 0,
        atol = 1e-6,
    )

    bus_slack_participation_factors = zeros(Float64, length(bus_numbers))
    for bn in bus_numbers
        bus = get_bus(sys, bn)
        idx = data.bus_lookup[bn]
        data.bus_type[idx, 1] == ACBusTypes.PQ && continue
        bus_slack_participation_factors[idx] = data.bus_slack_participation_factors[idx, 1]
    end

    # needed for the test implementation when data.bus_slack_participation_factors is compared to bus_slack_participation_factors
    generator_slack_participation_factors2 = Dict(
        (typeof(x), get_name(x)) => f for
        (x, f) in zip(devices, factors) if get_bustype(get_bus(x)) != ACBusTypes.PQ
    )
    pf2 = ACPowerFlow(ACSolver;
        generator_slack_participation_factors = generator_slack_participation_factors2,
    )

    _check_ds_pf(
        pf2,
        sys,
        bus_slack_participation_factors,
        bus_numbers,
        original_bus_power,
        original_gen_power,
        data_original_bus_power;
        check_connectivity = false,
    )
end

@testset "AC PF DS with several REF buses" for ACSolver in (
    NewtonRaphsonACPowerFlow,
    TrustRegionACPowerFlow,
)
    sys = System(100.0)

    n = 4

    buses_ref = [_add_simple_bus!(sys, i, ACBusTypes.REF, 230, 1.1, 0.0) for i in 1:n]
    buses_pv =
        [_add_simple_bus!(sys, i, ACBusTypes.PV, 230, 1.1, 0.0) for i in (n + 1):(2n)]
    buses_pq =
        [_add_simple_bus!(sys, i, ACBusTypes.PQ, 230, 1.1, 0.0) for i in (2n + 1):(3n)]

    loads = [_add_simple_load!(sys, b, 6, 2) for b in buses_pq]
    sources = [_add_simple_source!(sys, b, 0.0, 0.0) for b in buses_ref]
    gens = [_add_simple_thermal_standard!(sys, b, 0.0, 0.0) for b in buses_pv]
    all_buses = vcat(buses_ref, buses_pv, buses_pq)
    lines = [
        _add_simple_line!(sys, b1, b2, 1e-3, 1e-3, 0.0) for
        (b1, b2) in zip(all_buses[1:(end - 1)], all_buses[2:end])
    ]
    generator_slack_participation_factors =
        Dict((typeof(x), get_name(x)) => 1.0 for x in vcat(sources, gens))

    pf = ACPowerFlow(ACSolver;
        generator_slack_participation_factors = generator_slack_participation_factors,
    )
    solve_powerflow!(pf, sys; correct_bustypes = true)

    # equal slack participation
    for (s, g) in zip(sources, gens)
        @test isapprox(get_active_power(s), get_active_power(g), rtol = 0, atol = 1e-6)
    end

    for (g1, g2) in zip(gens[1:(end - 1)], gens[2:end])
        @test isapprox(get_active_power(g1), get_active_power(g2), rtol = 0, atol = 1e-6)
    end

    # test that isolated islands raise error
    b = _add_simple_bus!(sys, 100, ACBusTypes.PQ, 230, 1.1, 0.0)

    @test_throws "No REF bus found in the subnetwork" solve_powerflow!(
        pf,
        sys;
        check_connectivity = false,
        correct_bustypes = true,
    )
end

@testset "AC PF 10k bus system: voltage magnitudes" begin
    sys = PSB.build_system(
        PSB.MatpowerTestSystems,
        "matpower_ACTIVSg10k_sys";
        force_build = false,
    )
    @assert !isempty(get_components(PhaseShiftingTransformer, sys)) "System should have " *
                                                                    "phase shifting transformers: " *
                                                                    "change `force_build` to `true` in the test."
    pf_tr = ACPowerFlow{TrustRegionACPowerFlow}()
    data_tr = PowerFlowData(pf_tr, sys; correct_bustypes = true)
    solve_powerflow!(
        data_tr;
        pf = pf_tr,
        maxIterations = 200,
        factor = 0.1,
    )
    @test all(data_tr.bus_magnitude[:, 1] .<= 1.1)
    @test all(data_tr.bus_magnitude[:, 1] .>= 0.9)
end

function PowerFlowData_to_DataFrame(data::PowerFlowData)
    nbuses = size(data.bus_magnitude, 1)
    # Convert the PowerFlowData to a DataFrame
    bus_rev_lookup = fill(-1, nbuses)
    for (bus_no, row_no) in data.bus_lookup
        bus_rev_lookup[row_no] = bus_no
    end
    @assert !(-1 in bus_rev_lookup)
    df = DataFrame(;
        bus_number = bus_rev_lookup,
        Vm = data.bus_magnitude[:, 1],
        bus_type = data.bus_type[:, 1],
        angle = data.bus_angles[:, 1],
        generator_p = data.bus_activepower_injection[:, 1],
        generator_q = data.bus_reactivepower_injection[:, 1],
        load_p = data.bus_activepower_withdrawals[:, 1],
        load_q = data.bus_reactivepower_withdrawals[:, 1],
    )
    sort!(df, :bus_number)
    return df
end

@testset "ACTIVSg2000 matches matpower's solution" begin
    MATPOWER_CSV = joinpath(TEST_FILES_DIR, "test_data", "ACTIVSg2000_solved.csv")
    matpower_df = DataFrame(CSV.File(MATPOWER_CSV))

    sys_sienna = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
    pf_sienna = ACPowerFlow()
    data_sienna = PowerFlowData(pf_sienna, sys_sienna; correct_bustypes = true)
    solve_powerflow!(data_sienna; pf = pf_sienna, tol = 1e-11)
    sienna_df = PowerFlowData_to_DataFrame(data_sienna)
    @assert all(sienna_df[!, "bus_number"] .== matpower_df[!, "bus_number"])
    # The bus types don't match, so we don't compare them here. (We changed PQ with 
    # generators to PV. Matpower doesn't do this, though it treats them as PV internally.)
    @test norm(sienna_df[!, "Vm"] .- matpower_df[!, "Vm"], Inf) < 1e-3
    @test norm(sienna_df[!, "angle"] .- matpower_df[!, "angle"], Inf) < 1e-3
    @test norm(sienna_df[!, "generator_p"] .- matpower_df[!, "generator_p"], Inf) < 1e-3
    # this fails with 1e-4.
    @test norm(sienna_df[!, "generator_q"] .- matpower_df[!, "generator_q"], Inf) < 1e-3
end

#=
# Unfortunately, we correct the voltage magnitude to something sensible before running the 
# solver, so testing this isn't straightforward.
@testset "voltage validation" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    pf = ACPowerFlow{TrustRegionACPowerFlow}()
    data = PowerFlowData(pf, sys; correct_bustypes = true)
    solve_powerflow!(data; pf = pf)
    data.bus_magnitude[1, 1] = 2.0
    @test_logs (:warn, r".*voltage magnitudes outside of range.*") match_mode = :any solve_powerflow!(data; pf = pf)
end=#

@testset "Test ZIP loads: constant current" begin
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 2, ACBusTypes.PQ, 230, 1.1, 0.0)
    l = _add_simple_line!(sys, b1, b2, 5e-3, 5e-3, 1e-3)
    s1 = _add_simple_source!(sys, b1, 0.0, 0.0)
    lc = _add_simple_zip_load!(
        sys,
        b2;
        constant_current_active_power = 2.0,
        constant_current_reactive_power = 1.0,
    )
    data = PowerFlowData(ACPowerFlow(), sys; correct_bustypes = true)
    solve_powerflow!(data)

    # here we calculate the current that is observed in the power flow solution
    # supplied through the line to the load, which is a constant current load.
    s_t =
        data.branch_activepower_flow_to_from + 1im * data.branch_reactivepower_flow_to_from
    i_t = abs(s_t[1]) / data.bus_magnitude[2, 1] / sqrt(3)

    # get the load inputs from the load component
    load_input_power = (get_current_active_power(lc) + 1im * get_current_reactive_power(lc))
    # calculating by hand the current that corresponds to the load inputs
    # constant current load is given for 1.0 p.u. base voltage:
    load_input_current = abs(load_input_power) / 1.0 / sqrt(3)

    # Calculate the observed bus power based on the Ybus and voltage at the bus of the load:
    V = data.bus_magnitude[:, 1] .* exp.(1im * data.bus_angles[:, 1])
    Sbus = V .* conj(data.power_network_matrix.data * V)

    # The expected power at the bus, corresponding to the formula for constant current loads:
    load_expected_power = load_input_power * data.bus_magnitude[2, 1]

    # - due to the reference frame in Sbus vs. load inputs
    @test isapprox(Sbus[2], -load_expected_power; atol = 1e-6, rtol = 0)

    # Calculate the expected current based on the load inputs and the bus voltage:
    load_expected_current = abs(load_expected_power) / data.bus_magnitude[2, 1] / sqrt(3)

    @test isapprox(load_expected_current, load_input_current; atol = 1e-6, rtol = 0)

    @test isapprox(i_t, load_input_current; atol = 1e-6, rtol = 0)

    @test isapprox(data.bus_activepower_injection[2, 1], 0.0, atol = 1e-12, rtol = 0)
    @test isapprox(data.bus_reactivepower_injection[2, 1], 0.0, atol = 1e-12, rtol = 0)
end

@testset "Test ZIP loads: constant impedance" begin
    sys = System(100.0)
    b_1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 2, ACBusTypes.PQ, 230, 1.1, 0.0)
    l = _add_simple_line!(sys, b_1, b2, 5e-3, 5e-3, 1e-3)
    s1 = _add_simple_source!(sys, b_1, 0.0, 0.0)
    lz = _add_simple_zip_load!(
        sys,
        b2;
        constant_impedance_active_power = 2.0,
        constant_impedance_reactive_power = 1.0,
    )
    data = PowerFlowData(ACPowerFlow(), sys; correct_bustypes = true)
    solve_powerflow!(data)

    # here we calculate the current that is observed in the power flow solution
    # supplied through the line to the load, which is a constant impedance load.
    s_t =
        data.branch_activepower_flow_to_from + 1im * data.branch_reactivepower_flow_to_from
    i_t = abs(s_t[1]) / data.bus_magnitude[2, 1] / sqrt(3)

    # get the load inputs from the load component
    load_input_power =
        (get_impedance_active_power(lz) + 1im * get_impedance_reactive_power(lz))
    # calculating by hand the impedance that corresponds to the load inputs
    # constant impedance load is given for 1.0 p.u. base voltage:
    load_input_impedance = 1.0^2 / abs(load_input_power)

    # Calculate the observed bus power based on the Ybus and voltage at the bus of the load:
    V = data.bus_magnitude[:, 1] .* exp.(1im * data.bus_angles[:, 1])
    Sbus = V .* conj(data.power_network_matrix.data * V)

    # The expected power at the bus, corresponding to the formula for constant impedance loads:
    load_expected_power = load_input_power * data.bus_magnitude[2, 1]^2

    # - due to the reference frame in Sbus vs. load inputs
    @test isapprox(Sbus[2], -load_expected_power; atol = 1e-6, rtol = 0)

    # Calculate the expected impedance based on the load inputs and the bus voltage:
    load_expected_impedance = (data.bus_magnitude[2, 1]^2) / abs(load_expected_power)

    @test isapprox(load_expected_impedance, load_input_impedance; atol = 1e-6, rtol = 0)

    s_zip_load =
        PF.get_bus_activepower_total_withdrawals(data, 2, 1) +
        1im * PF.get_bus_reactivepower_total_withdrawals(data, 2, 1)

    @test isapprox(
        Sbus[2],
        -s_zip_load;
        atol = 1e-6,
        rtol = 0,
    )

    @test isapprox(
        data.bus_magnitude[2, 1]^2 / load_input_impedance,
        abs(s_zip_load);
        atol = 1e-6,
        rtol = 0,
    )

    @test isapprox(  # <- need to come back to this and check why it fails
        s_t[1],
        -s_zip_load;
        atol = 1e-6,
        rtol = 0,
    )

    @test isapprox(data.bus_activepower_injection[2, 1], 0.0, atol = 1e-12, rtol = 0)
    @test isapprox(data.bus_reactivepower_injection[2, 1], 0.0, atol = 1e-12, rtol = 0)
end
