@testset "TimePowerFlowData" begin
    n_buses = 5
    n_arcs = 4
    n_time_steps = 3

    @testset "Default construction" begin
        results = TimePowerFlowData(n_buses, n_arcs, n_time_steps)

        @test size(results.bus_magnitude) == (n_buses, n_time_steps)
        @test size(results.bus_angles) == (n_buses, n_time_steps)
        @test size(results.bus_type) == (n_buses, n_time_steps)
        @test size(results.bus_active_power_injections) == (n_buses, n_time_steps)
        @test size(results.arc_active_power_flow_from_to) == (n_arcs, n_time_steps)
        @test length(results.converged) == n_time_steps

        @test results.loss_factors === nothing
        @test results.voltage_stability_factors === nothing
        @test results.arc_active_power_losses === nothing
        @test results.lcc_active_power_flow_from_to === nothing
        @test results.lcc_active_power_flow_to_from === nothing
        @test results.lcc_reactive_power_flow_from_to === nothing
        @test results.lcc_reactive_power_flow_to_from === nothing
    end

    @testset "Optional fields enabled" begin
        results = TimePowerFlowData(
            n_buses,
            n_arcs,
            n_time_steps;
            calculate_loss_factors = true,
            calculate_voltage_stability_factors = true,
            make_arc_active_power_losses = true,
        )

        @test size(results.loss_factors) == (n_buses, n_time_steps)
        @test size(results.voltage_stability_factors) == (n_buses, n_time_steps)
        @test size(results.arc_active_power_losses) == (n_arcs, n_time_steps)
    end

    @testset "LCC fields" begin
        n_lccs = 2
        results = TimePowerFlowData(n_buses, n_arcs, n_time_steps; n_lccs = n_lccs)

        @test size(results.lcc_active_power_flow_from_to) == (n_lccs, n_time_steps)
        @test size(results.lcc_active_power_flow_to_from) == (n_lccs, n_time_steps)
        @test size(results.lcc_reactive_power_flow_from_to) == (n_lccs, n_time_steps)
        @test size(results.lcc_reactive_power_flow_to_from) == (n_lccs, n_time_steps)
    end
end

@testset "TimeContingencyPowerFlowData" begin
    n_buses = 5
    n_arcs = 4
    n_time_steps = 3
    ctg_labels = ["base", "line_1_out", "line_2_out"]

    @testset "Default construction" begin
        results = TimeContingencyPowerFlowData(n_buses, n_arcs, n_time_steps, ctg_labels)
        n_ctg = length(ctg_labels)

        @test size(results.bus_magnitude) == (n_buses, n_time_steps, n_ctg)
        @test size(results.bus_active_power_injections) == (n_buses, n_time_steps, n_ctg)
        @test size(results.arc_active_power_flow_from_to) == (n_arcs, n_time_steps, n_ctg)
        @test size(results.converged) == (n_time_steps, n_ctg)

        @test results.loss_factors === nothing
        @test results.voltage_stability_factors === nothing
        @test results.arc_active_power_losses === nothing
        @test results.lcc_active_power_flow_from_to === nothing
        @test results.lcc_active_power_flow_to_from === nothing
        @test results.lcc_reactive_power_flow_from_to === nothing
        @test results.lcc_reactive_power_flow_to_from === nothing
        @test all(isnothing, results.network_modifications)
    end

    @testset "Contingency lookup" begin
        results = TimeContingencyPowerFlowData(n_buses, n_arcs, n_time_steps, ctg_labels)
        lookup = PowerFlows.get_contingency_lookup(results)

        @test lookup["base"] == 1
        @test lookup["line_1_out"] == 2
        @test lookup["line_2_out"] == 3
        @test PowerFlows.get_contingency_labels(results) == ctg_labels
        @test PowerFlows.get_n_contingencies(results) == 3
    end

    @testset "NetworkModification storage and access" begin
        mod1 = PNM.NetworkModification(
            "line_1_out",
            [PNM.ArcModification(1, -0.5)],
        )
        mod2 = PNM.NetworkModification(
            "line_2_out",
            [PNM.ArcModification(2, -0.3)],
        )
        mods = Union{Nothing, PNM.NetworkModification}[nothing, mod1, mod2]
        results = TimeContingencyPowerFlowData(
            n_buses,
            n_arcs,
            n_time_steps,
            ctg_labels;
            network_modifications = mods,
        )

        @test PowerFlows.get_network_modification(results, "base") === nothing
        @test PowerFlows.get_network_modification(results, "line_1_out") === mod1
        @test PowerFlows.get_network_modification(results, "line_2_out") === mod2
        @test PowerFlows.get_network_modification(results, 1) === nothing
        @test PowerFlows.get_network_modification(results, 2) === mod1
        @test length(PowerFlows.get_network_modifications(results)) == 3
    end

    @testset "set_network_modification!" begin
        results = TimeContingencyPowerFlowData(n_buses, n_arcs, n_time_steps, ctg_labels)
        @test PowerFlows.get_network_modification(results, "base") === nothing

        new_mod = PNM.NetworkModification(
            "base_mod",
            [PNM.ArcModification(3, -0.1)],
        )
        PowerFlows.set_network_modification!(results, "base", new_mod)
        @test PowerFlows.get_network_modification(results, "base") === new_mod
    end

    @testset "Optional fields enabled" begin
        results = TimeContingencyPowerFlowData(
            n_buses,
            n_arcs,
            n_time_steps,
            ctg_labels;
            calculate_loss_factors = true,
            calculate_voltage_stability_factors = true,
            make_arc_active_power_losses = true,
        )
        n_ctg = length(ctg_labels)

        @test size(results.loss_factors) == (n_buses, n_time_steps, n_ctg)
        @test size(results.voltage_stability_factors) == (n_buses, n_time_steps, n_ctg)
        @test size(results.arc_active_power_losses) == (n_arcs, n_time_steps, n_ctg)
    end

    @testset "LCC fields" begin
        n_lccs = 2
        results = TimeContingencyPowerFlowData(
            n_buses, n_arcs, n_time_steps, ctg_labels; n_lccs = n_lccs,
        )
        n_ctg = length(ctg_labels)

        @test size(results.lcc_active_power_flow_from_to) == (n_lccs, n_time_steps, n_ctg)
        @test size(results.lcc_active_power_flow_to_from) == (n_lccs, n_time_steps, n_ctg)
        @test size(results.lcc_reactive_power_flow_from_to) ==
              (n_lccs, n_time_steps, n_ctg)
        @test size(results.lcc_reactive_power_flow_to_from) ==
              (n_lccs, n_time_steps, n_ctg)
    end

    @testset "Constructor validation" begin
        bad_mods = Union{Nothing, PNM.NetworkModification}[nothing, nothing]
        @test_throws ArgumentError TimeContingencyPowerFlowData(
            n_buses, n_arcs, n_time_steps, ctg_labels;
            network_modifications = bad_mods,
        )
        @test_throws ArgumentError TimeContingencyPowerFlowData(
            n_buses, n_arcs, n_time_steps, ["base", "ctg_1", "base"],
        )
    end
end

@testset "PowerFlowData property forwarding" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    pf = ACPowerFlow()
    data = PowerFlows.PowerFlowData(pf, sys)

    @test data.results isa PowerFlows.TimePowerFlowData
    @test data.bus_magnitude === data.results.bus_magnitude
    @test data.bus_active_power_injections === data.results.bus_active_power_injections
    @test data.converged === data.results.converged

    r = PowerFlows.get_results(data)
    @test r === data.results
end

@testset "Contingency-aware slicing" begin
    n_buses = 4
    n_arcs = 5
    n_time_steps = 2
    labels = ["base", "ctg_1", "ctg_2"]

    results = PowerFlows.TimeContingencyPowerFlowData(
        n_buses, n_arcs, n_time_steps, labels,
    )
    results.bus_magnitude[3, 1, 2] = 0.97

    slice = PowerFlows.get_contingency_slice(results, :bus_magnitude, 2)
    @test size(slice) == (n_buses, n_time_steps)
    @test slice[3, 1] == 0.97

    slice2 = PowerFlows.get_contingency_slice(results, :bus_magnitude, "ctg_1")
    @test slice2 == slice

    # Verify it's a view (mutation reflects back)
    slice[1, 1] = 0.5
    @test results.bus_magnitude[1, 1, 2] == 0.5
end

@testset "get_contingency_slice validation" begin
    results = PowerFlows.TimeContingencyPowerFlowData(
        4, 5, 2, ["base", "ctg_1", "ctg_2"],
    )

    @test_throws ArgumentError PowerFlows.get_contingency_slice(
        results, :nonexistent_field, 1)
    @test_throws ArgumentError PowerFlows.get_contingency_slice(
        results, :loss_factors, 1)
    @test_throws ArgumentError PowerFlows.get_contingency_slice(
        results, :bus_magnitude, 0)
    @test_throws ArgumentError PowerFlows.get_contingency_slice(
        results, :bus_magnitude, 4)
    @test_throws ArgumentError PowerFlows.get_contingency_slice(
        results, :bus_magnitude, "nonexistent")
end
