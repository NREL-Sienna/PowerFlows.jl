@testset "TimePowerFlowData" begin
    n_buses = 5
    n_arcs = 4
    n_time_steps = 3

    @testset "Default construction" begin
        results = TimePowerFlowData(n_buses, n_arcs, n_time_steps)

        # Bus fields have correct dimensions.
        @test size(results.bus_magnitude) == (n_buses, n_time_steps)
        @test size(results.bus_angles) == (n_buses, n_time_steps)
        @test size(results.bus_type) == (n_buses, n_time_steps)

        # Bus magnitude defaults to ones (flat start).
        @test all(results.bus_magnitude .== 1.0)
        # Bus angles default to zeros.
        @test all(results.bus_angles .== 0.0)
        # Bus types default to PQ.
        @test all(results.bus_type .== Ref(PowerSystems.ACBusTypes.PQ))

        # Arc fields have correct dimensions and default to zeros.
        for field in [
            :arc_active_power_flow_from_to,
            :arc_reactive_power_flow_from_to,
            :arc_active_power_flow_to_from,
            :arc_reactive_power_flow_to_from,
            :arc_angle_differences,
        ]
            mat = getfield(results, field)
            @test size(mat) == (n_arcs, n_time_steps)
            @test all(mat .== 0.0)
        end

        # Converged defaults to all false.
        @test length(results.converged) == n_time_steps
        @test !any(results.converged)

        # Optional fields default to nothing.
        @test results.loss_factors === nothing
        @test results.voltage_stability_factors === nothing
        @test results.arc_active_power_losses === nothing
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

        @test results.loss_factors !== nothing
        @test size(results.loss_factors) == (n_buses, n_time_steps)
        @test all(results.loss_factors .== 0.0)

        @test results.voltage_stability_factors !== nothing
        @test size(results.voltage_stability_factors) == (n_buses, n_time_steps)
        @test all(results.voltage_stability_factors .== 0.0)

        @test results.arc_active_power_losses !== nothing
        @test size(results.arc_active_power_losses) == (n_arcs, n_time_steps)
        @test all(results.arc_active_power_losses .== 0.0)
    end

    @testset "Subtype relationship" begin
        @test TimePowerFlowData <: AbstractPowerFlowResults
    end

    @testset "Single time step" begin
        results = TimePowerFlowData(10, 8, 1)
        @test size(results.bus_magnitude) == (10, 1)
        @test size(results.arc_active_power_flow_from_to) == (8, 1)
        @test length(results.converged) == 1
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

        # Bus fields have correct 3D dimensions (entity, time_step, contingency).
        @test size(results.bus_magnitude) == (n_buses, n_time_steps, n_ctg)
        @test size(results.bus_angles) == (n_buses, n_time_steps, n_ctg)
        @test size(results.bus_type) == (n_buses, n_time_steps, n_ctg)

        # Bus magnitude defaults to ones (flat start).
        @test all(results.bus_magnitude .== 1.0)
        # Bus angles default to zeros.
        @test all(results.bus_angles .== 0.0)
        # Bus types default to PQ.
        @test all(results.bus_type .== Ref(PowerSystems.ACBusTypes.PQ))

        # Arc fields have correct 3D dimensions and default to zeros.
        for field in [
            :arc_active_power_flow_from_to,
            :arc_reactive_power_flow_from_to,
            :arc_active_power_flow_to_from,
            :arc_reactive_power_flow_to_from,
            :arc_angle_differences,
        ]
            arr = getfield(results, field)
            @test size(arr) == (n_arcs, n_time_steps, n_ctg)
            @test all(arr .== 0.0)
        end

        # Converged is a BitMatrix of (time_steps, contingencies), all false.
        @test size(results.converged) == (n_time_steps, n_ctg)
        @test !any(results.converged)

        # Optional fields default to nothing.
        @test results.loss_factors === nothing
        @test results.voltage_stability_factors === nothing
        @test results.arc_active_power_losses === nothing

        # Network modifications default to nothing for each contingency.
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

        # Access by label.
        @test PowerFlows.get_network_modification(results, "base") === nothing
        @test PowerFlows.get_network_modification(results, "line_1_out") === mod1
        @test PowerFlows.get_network_modification(results, "line_2_out") === mod2

        # Access by index.
        @test PowerFlows.get_network_modification(results, 1) === nothing
        @test PowerFlows.get_network_modification(results, 2) === mod1

        # get_network_modifications returns the full vector.
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

        @test results.loss_factors !== nothing
        @test size(results.loss_factors) == (n_buses, n_time_steps, n_ctg)
        @test all(results.loss_factors .== 0.0)

        @test results.voltage_stability_factors !== nothing
        @test size(results.voltage_stability_factors) == (n_buses, n_time_steps, n_ctg)
        @test all(results.voltage_stability_factors .== 0.0)

        @test results.arc_active_power_losses !== nothing
        @test size(results.arc_active_power_losses) == (n_arcs, n_time_steps, n_ctg)
        @test all(results.arc_active_power_losses .== 0.0)
    end

    @testset "Subtype relationship" begin
        @test TimeContingencyPowerFlowData <: AbstractPowerFlowResults
    end

    @testset "Mismatched network_modifications length throws" begin
        bad_mods = Union{Nothing, PNM.NetworkModification}[nothing, nothing]
        @test_throws ArgumentError TimeContingencyPowerFlowData(
            n_buses,
            n_arcs,
            n_time_steps,
            ctg_labels;
            network_modifications = bad_mods,
        )
    end

    @testset "Single contingency" begin
        results =
            TimeContingencyPowerFlowData(n_buses, n_arcs, n_time_steps, ["base_case"])
        @test size(results.bus_magnitude) == (n_buses, n_time_steps, 1)
        @test size(results.converged) == (n_time_steps, 1)
        @test PowerFlows.get_n_contingencies(results) == 1
    end
end

@testset "PowerFlowData contains TimePowerFlowData results" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    pf = ACPowerFlow()
    data = PowerFlows.PowerFlowData(pf, sys)

    @test data.results isa PowerFlows.TimePowerFlowData
    @test data.bus_magnitude === data.results.bus_magnitude
    @test data.bus_angles === data.results.bus_angles
    @test data.bus_type === data.results.bus_type
    @test data.converged === data.results.converged

    # Mutation through forwarding works
    data.bus_magnitude[1, 1] = 0.95
    @test data.results.bus_magnitude[1, 1] == 0.95
end

@testset "get_results accessor" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    pf = ACPowerFlow()
    data = PowerFlows.PowerFlowData(pf, sys)

    r = PowerFlows.get_results(data)
    @test r isa PowerFlows.TimePowerFlowData
    @test r.bus_magnitude === data.bus_magnitude
end

@testset "Contingency-aware slicing on TimeContingencyPowerFlowData" begin
    n_buses = 4
    n_arcs = 5
    n_time_steps = 2
    labels = ["base", "ctg_1", "ctg_2"]

    results = PowerFlows.TimeContingencyPowerFlowData(
        n_buses, n_arcs, n_time_steps, labels,
    )
    # Write a value into contingency 2, time_step 1, bus 3
    results.bus_magnitude[3, 1, 2] = 0.97

    # Slice by contingency index
    slice = PowerFlows.get_contingency_slice(results, :bus_magnitude, 2)
    @test size(slice) == (n_buses, n_time_steps)
    @test slice[3, 1] == 0.97

    # Slice by contingency label
    slice2 = PowerFlows.get_contingency_slice(results, :bus_magnitude, "ctg_1")
    @test slice2 === slice

    # Verify it's a view (mutation reflects back)
    slice[1, 1] = 0.5
    @test results.bus_magnitude[1, 1, 2] == 0.5
end
