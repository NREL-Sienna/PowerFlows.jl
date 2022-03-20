result = [
    2.32551,
    -0.155293,
    0.469214,
    -0.0870457,
    0.271364,
    -0.222398,
    1.01423,
    -0.179009,
    1.01724,
    -0.152972,
    0.216039,
    -0.251637,
    1.05034,
    -0.231289,
    0.245388,
    -0.231289,
    1.03371,
    -0.258872,
    1.03256,
    -0.262519,
    1.04748,
    -0.259143,
    1.0535,
    -0.266484,
    1.04711,
    -0.267177,
    1.02131,
    -0.280381,
]

p_gen_matpower_3bus = [20.3512373930753, 100.0, 100.0]
q_gen_matpower_3bus = [45.516916781567232, 10.453799727283879, -31.992561631394636]

pf_sys5_re = PSB.build_system(PSB.PSITestSystems, "c_sys5_re"; add_forecasts=false)
remove_component!(Line, pf_sys5_re, "1")
remove_component!(Line, pf_sys5_re, "2")
br = get_component(Line, pf_sys5_re, "6")
PSY.set_x!(br, 20.0)
PSY.set_r!(br, 2.0)

@testset "NLsolve Power Flow testing" begin
    # This is a negative test. The data passed for sys5_re is known to be infeasible.
    @test_logs(
        (:error, "The powerflow solver returned convergence = false"),
        match_mode = :any,
        @test !run_powerflow!(pf_sys5_re, finite_diff=true)
    )
    #Compare results between finite diff methods and Jacobian method
    res_finite_diff = run_powerflow(
        PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts=false),
        finite_diff=true,
    )
    res_jacobian =
        run_powerflow(PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts=false))
    @test LinearAlgebra.norm(
        res_finite_diff["bus_results"].Vm - res_jacobian["bus_results"].Vm,
    ) <= 1e-6
    @test run_powerflow!(
        PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts=false),
        finite_diff=true,
        method=:newton,
    )

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts=false)
    branch = first(PSY.get_components(Line, sys))
    dyn_branch = DynamicBranch(branch)
    add_component!(sys, dyn_branch)
    @test dyn_pf = run_powerflow!(sys)
    dyn_pf = run_powerflow(sys)
    @test LinearAlgebra.norm(dyn_pf["bus_results"].Vm - res_jacobian["bus_results"].Vm) <=
          1e-6

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts=false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    run_powerflow!(sys)
    @test PSY.get_active_power_flow(line) == 0.0
    test_bus = get_component(PSY.Bus, sys, "Bus 4")
    @test isapprox(PSY.get_magnitude(test_bus), 1.002; atol=1e-3)

    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts=false)
    line = get_component(Line, sys, "Line4")
    PSY.set_available!(line, false)
    res = run_powerflow(sys)
    @test res["flow_results"].P_from_to[4] == 0.0
    @test res["flow_results"].P_to_from[4] == 0.0

    sys = PSB.build_system(PSB.PSSETestSystems, "psse_240_case_renewable_sys")
    @test run_powerflow!(sys)

    sys_3bus = PSB.build_system(PSB.PSSETestSystems, "psse_3bus_gen_cls_sys")
    bus_103 = get_component(PSY.Bus, sys_3bus, "BUS 3")
    fix_shunt = PSY.FixedAdmittance("FixAdmBus3", true, bus_103, 0.0 + 0.2im)
    add_component!(sys_3bus, fix_shunt)
    df = run_powerflow(sys_3bus)
    @test isapprox(df["bus_results"].P_gen, p_gen_matpower_3bus, atol=1e-4)
    @test isapprox(df["bus_results"].Q_gen, q_gen_matpower_3bus, atol=1e-4)
end

@testset "Test 240 Case PSS/e results" begin
    file = joinpath(
        TEST_FILES_DIR,
        "test_data",
        "WECC240_v04_DPV_RE20_v33_6302_xfmr_DPbuscode_PFadjusted_V32_noRemoteVctrl.raw",
    )
    system = System(
        file,
        bus_name_formatter=x -> strip(string(x["name"])) * "-" * string(x["index"]),
        runchecks=false,
    )

    pf_bus_result_file = joinpath(TEST_FILES_DIR, "test_data", "pf_bus_results.csv")
    pf_gen_result_file = joinpath(TEST_FILES_DIR, "test_data", "pf_gen_results.csv")

    pf = run_powerflow!(system)
    @test pf
    pf_result_df = run_powerflow(system)

    v_diff, angle_diff, number = psse_bus_results_compare(pf_bus_result_file, pf_result_df)
    p_diff, q_diff, names = psse_gen_results_compare(pf_gen_result_file, system)

    base_power = get_base_power(system)
    @test norm(v_diff, Inf) < DIFF_INF_TOLERANCE
    @test norm(v_diff, 2) / length(v_diff) < DIFF_L2_TOLERANCE
    @test norm(angle_diff, Inf) < DIFF_INF_TOLERANCE
    @test norm(angle_diff, 2) / length(angle_diff) < DIFF_L2_TOLERANCE
    @test norm(p_diff, Inf) < DIFF_INF_TOLERANCE*base_power
    @test norm(p_diff, 2) / length(p_diff) < DIFF_L2_TOLERANCE
    @test sum(q_diff) < DIFF_INF_TOLERANCE*base_power
    @test norm(q_diff, 2) / length(q_diff) < DIFF_L2_TOLERANCE
end
