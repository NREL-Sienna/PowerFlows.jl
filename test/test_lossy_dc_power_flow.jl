@testset "Lossy DCLF: default matches explicit lossy_flows=false" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    default_result = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    explicit_false = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true, lossy_flows = false),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    @test default_result["1"]["flow_results"].P_from_to ==
          explicit_false["1"]["flow_results"].P_from_to
end

@testset "Lossy DCLF: lossy_flows=true gives different and asymmetric flows" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    lossless = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    lossy = solve_power_flow(
        DCPowerFlow(; correct_bustypes = true, lossy_flows = true),
        sys,
        PF.FlowReporting.ARC_FLOWS,
    )
    # Flows should differ (c_sys14 has nonzero resistance).
    @test lossless["1"]["flow_results"].P_from_to != lossy["1"]["flow_results"].P_from_to
    # Lossy model reports asymmetric flows: P_to_from != -P_from_to
    lossy_flows = lossy["1"]["flow_results"]
    @test lossy_flows.P_to_from != -lossy_flows.P_from_to
    # Loss conservation: P_from_to + P_to_from ≈ P_losses
    @test isapprox(
        lossy_flows.P_from_to .+ lossy_flows.P_to_from,
        lossy_flows.P_losses;
        atol = 1e-6,
    )
    # Lossless model still has symmetric flows.
    lossless_flows = lossless["1"]["flow_results"]
    @test lossless_flows.P_to_from == -lossless_flows.P_from_to
end

@testset "Lossy DCLF: arc_lossy_admittance fields" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # Without lossy_flows: arc admittance fields should be nothing.
    data_no_loss = PowerFlowData(DCPowerFlow(; correct_bustypes = true), sys)
    @test data_no_loss.arc_lossy_admittance_from_to === nothing
    @test data_no_loss.arc_lossy_admittance_to_from === nothing

    # With lossy_flows: arc admittance fields should be set.
    data_lossy =
        PowerFlowData(DCPowerFlow(; correct_bustypes = true, lossy_flows = true), sys)
    @test data_lossy.arc_lossy_admittance_from_to !== nothing
    @test data_lossy.arc_lossy_admittance_to_from !== nothing
end

@testset "Lossy DCLF: multi-period" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    data = PowerFlowData(
        DCPowerFlow(;
            lossy_flows = true,
            time_steps = 3,
        ),
        sys,
    )
    # Replicate time step 1 data (injections/withdrawals) to all time steps.
    for col in 2:3
        data.bus_active_power_injections[:, col] .=
            data.bus_active_power_injections[:, 1]
        data.bus_active_power_withdrawals[:, col] .=
            data.bus_active_power_withdrawals[:, 1]
    end
    solve_power_flow!(data)
    # All time steps should produce identical results (same system state).
    @test data.arc_active_power_flow_from_to[:, 1] ==
          data.arc_active_power_flow_from_to[:, 2]
    @test data.arc_active_power_flow_from_to[:, 1] ==
          data.arc_active_power_flow_from_to[:, 3]
end

@testset "Lossy DCLF: AC and PTDF data have no lossy admittances" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    ac_data = PowerFlowData(ACPowerFlow(; correct_bustypes = true), sys)
    @test ac_data.arc_lossy_admittance_from_to === nothing
    ptdf_data = PowerFlowData(PTDFDCPowerFlow(; correct_bustypes = true), sys)
    @test ptdf_data.arc_lossy_admittance_from_to === nothing
end

function _build_arc_flow_map(data::PowerFlowData, sys::PSY.System)
    arc_tuples = PNM.get_arc_axis(data.aux_network_matrix)
    p_from_to = data.arc_active_power_flow_from_to[:, 1] .* PSY.get_base_power(sys)
    p_to_from = data.arc_active_power_flow_to_from[:, 1] .* PSY.get_base_power(sys)
    p_losses = data.arc_active_power_losses[:, 1] .* PSY.get_base_power(sys)
    return Dict(
        arc => (
            P_from_to = p_from_to[ix],
            P_to_from = p_to_from[ix],
            P_losses = p_losses[ix],
        ) for (ix, arc) in enumerate(arc_tuples)
    )
end

function _compare_arc_flows_to_psse_csv(
    data::PowerFlowData,
    sys::PSY.System,
    psse_csv_path::String;
    atol::Float64 = 1e-3,
)
    @test isfile(psse_csv_path)
    psse_df = CSV.read(psse_csv_path, DataFrame)
    arc_map = _build_arc_flow_map(data, sys)

    @test nrow(psse_df) == length(arc_map)

    for row in eachrow(psse_df)
        from_bus = row[Symbol("FromBus#")]
        to_bus = row[Symbol("ToBus#")]
        p_ij = row[Symbol("Pij(MW)")]
        p_ji = row[Symbol("Pji(MW)")]
        p_loss = row[Symbol("Ploss(MW)")]

        if haskey(arc_map, (from_bus, to_bus))
            arc = arc_map[(from_bus, to_bus)]
            @test isapprox(arc.P_from_to, p_ij; atol = atol)
            @test isapprox(arc.P_to_from, p_ji; atol = atol)
            @test isapprox(arc.P_losses, p_loss; atol = atol)
        elseif haskey(arc_map, (to_bus, from_bus))
            # If the arc orientation differs, compare against swapped Pij/Pji.
            arc = arc_map[(to_bus, from_bus)]
            @test isapprox(arc.P_from_to, p_ji; atol = atol)
            @test isapprox(arc.P_to_from, p_ij; atol = atol)
            @test isapprox(arc.P_losses, p_loss; atol = atol)
        else
            @test false
        end
    end
end

function _get_custom_uc_template_simulation!(template)
    PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalBasicDispatch)
    PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
    PSI.set_device_model!(template, PSY.Line, PSI.StaticBranchBounds)
    return
end

function _run_custom_pcm(sys::PSY.System, s_time::Dates.DateTime)
    PSY.transform_single_time_series!(sys, Dates.Hour(1), Dates.Hour(1))
    ptdf = PNM.PTDF(sys)

    power_flow_model = DCPowerFlow(; lossy_flows = true)
    template = PSI.ProblemTemplate(
        PSI.NetworkModel(
            PSI.PTDFPowerModel;
            use_slacks = false,
            power_flow_evaluation = power_flow_model,
            PTDF_matrix = ptdf,
        ),
    )

    _get_custom_uc_template_simulation!(template)

    highs_optimizer = JuMP.optimizer_with_attributes(
        HiGHS.Optimizer,
        "time_limit" => 100.0,
        "random_seed" => 12345,
        "log_to_console" => false,
    )

    return PSI.DecisionModel(
        template,
        sys;
        name = "UC1",
        optimizer = highs_optimizer,
        store_variable_names = true,
        initialize_model = false,
        optimizer_solve_log_print = true,
        direct_mode_optimizer = true,
        check_numerical_bounds = false,
        rebuild_model = false,
        system_to_file = true,
        initial_time = s_time,
    )
end

function _get_pcm_power_flow_data(model)
    container = PSI.get_optimization_container(model)
    pf_eval_data = PSI.get_power_flow_evaluation_data(container)
    @test !isempty(pf_eval_data)
    return PSI.get_power_flow_data(pf_eval_data[1])
end

@testset "Lossy DCLF: PCM arc flows match PSS/e results" begin
    PSB.clear_all_serialized_systems()
    pscb_case_name = "2Area 5 Bus System"
    sys = PSB.build_system(PSB.PSISystems, pscb_case_name)
    period = Dates.DateTime("2020-01-01T00:00:00")

    model = _run_custom_pcm(sys, period)
    PSI.build!(model; console_level = Logging.Info, output_dir = mktempdir())
    PSI.solve!(model)

    data = _get_pcm_power_flow_data(model)

    psse_csv_path = joinpath(
        TEST_DATA_DIR,
        "2area-5-bus-system",
        "export_1_1_psse_dcpf_results.csv",
    )
    _compare_arc_flows_to_psse_csv(data, sys, psse_csv_path)
end
