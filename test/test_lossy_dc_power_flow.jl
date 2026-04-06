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

function _compute_net_injections_from_psse_csv(
    psse_csv_path::String,
    bus_lookup::Dict{Int, Int},
    base_power::Float64,
)
    psse_df = CSV.read(psse_csv_path, DataFrame)
    n_buses = length(bus_lookup)
    net_injections = zeros(Float64, n_buses)
    for row in eachrow(psse_df)
        from_bus = row[Symbol("FromBus#")]
        to_bus = row[Symbol("ToBus#")]
        p_ij = row[Symbol("Pij(MW)")]
        p_ji = row[Symbol("Pji(MW)")]
        if haskey(bus_lookup, from_bus)
            net_injections[bus_lookup[from_bus]] += p_ij
        end
        if haskey(bus_lookup, to_bus)
            net_injections[bus_lookup[to_bus]] += p_ji
        end
    end
    net_injections ./= base_power
    return net_injections
end

@testset "Lossy DCLF: 2-area system arc flows match PSS/e results" begin
    sys = PSB.build_system(PSB.PSISystems, "2Area 5 Bus System")
    data = PowerFlowData(
        DCPowerFlow(; correct_bustypes = true, lossy_flows = true),
        sys,
    )

    psse_csv_path = joinpath(
        TEST_DATA_DIR,
        "2area-5-bus-system",
        "export_1_1_psse_dcpf_results.csv",
    )
    bus_lookup = PNM.get_bus_lookup(data.aux_network_matrix)
    net_injections =
        _compute_net_injections_from_psse_csv(
            psse_csv_path,
            bus_lookup,
            PSY.get_base_power(sys),
        )
    data.bus_active_power_injections[:, 1] .= max.(net_injections, 0.0)
    data.bus_active_power_withdrawals[:, 1] .= max.(-net_injections, 0.0)

    solve_power_flow!(data)
    _compare_arc_flows_to_psse_csv(data, sys, psse_csv_path; atol = 5.0)
end
