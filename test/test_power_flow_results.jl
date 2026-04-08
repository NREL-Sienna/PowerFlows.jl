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
