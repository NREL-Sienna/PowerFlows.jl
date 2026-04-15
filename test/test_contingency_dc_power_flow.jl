@testset "Contingency DC PF — N-1 angles match rebuild truth on c_sys5" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    branches = collect(
        get_components(
            x -> !(typeof(x) <: Union{PhaseShiftingTransformer, DiscreteControlledACBranch}),
            ACTransmission,
            sys,
        ),
    )
    outages = PSY.Outage[]
    for br in branches
        o = GeometricDistributionForcedOutage(;
            mean_time_to_recovery = 0.0,
            outage_transition_probability = 0.0,
        )
        add_supplemental_attribute!(sys, br, o)
        push!(outages, o)
    end

    result = solve_contingency_dc_power_flow(sys, outages; time_steps = 1)

    @test length(PowerFlows.get_contingency_labels(result)) == length(branches) + 1
    @test all(result.converged)

    for (k, br) in enumerate(branches)
        set_available!(br, false)
        ref_data = PowerFlowData(DCPowerFlow(; time_steps = 1), sys)
        solve_power_flow!(ref_data)
        set_available!(br, true)

        @test isapprox(
            result.bus_angles[:, :, k + 1],
            ref_data.bus_angles;
            atol = 1e-8,
        )
    end
end

@testset "Contingency DC PF — outaged arc flow is zero" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    line1 = get_component(Line, sys, "1")
    outage = GeometricDistributionForcedOutage(;
        mean_time_to_recovery = 0.0,
        outage_transition_probability = 0.0,
    )
    add_supplemental_attribute!(sys, line1, outage)

    result = solve_contingency_dc_power_flow(sys, [outage]; time_steps = 1)
    mod = PowerFlows.get_network_modification(result, 2)
    @test mod !== nothing
    for arc_mod in mod.arc_modifications
        @test isapprox(
            result.arc_active_power_flow_from_to[arc_mod.arc_index, 1, 2],
            0.0;
            atol = 1e-12,
        )
        @test isapprox(
            result.arc_active_power_flow_to_from[arc_mod.arc_index, 1, 2],
            0.0;
            atol = 1e-12,
        )
    end
end

@testset "Contingency DC PF — N-2 via shared Outage matches rebuild" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    line1 = get_component(Line, sys, "1")
    line2 = get_component(Line, sys, "2")
    outage = GeometricDistributionForcedOutage(;
        mean_time_to_recovery = 0.0,
        outage_transition_probability = 0.0,
    )
    add_supplemental_attribute!(sys, line1, outage)
    add_supplemental_attribute!(sys, line2, outage)

    result = solve_contingency_dc_power_flow(sys, [outage]; time_steps = 1)
    @test all(result.converged)

    set_available!(line1, false)
    set_available!(line2, false)
    ref_data = PowerFlowData(DCPowerFlow(; time_steps = 1), sys)
    solve_power_flow!(ref_data)
    set_available!(line1, true)
    set_available!(line2, true)

    @test isapprox(result.bus_angles[:, :, 2], ref_data.bus_angles; atol = 1e-8)
end

@testset "Contingency DC PF — Woodbury (:auto) matches :refactor" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    branches = collect(
        get_components(
            x -> !(typeof(x) <: Union{PhaseShiftingTransformer, DiscreteControlledACBranch}),
            ACTransmission,
            sys,
        ),
    )
    outages = PSY.Outage[]
    for br in branches
        o = GeometricDistributionForcedOutage(;
            mean_time_to_recovery = 0.0,
            outage_transition_probability = 0.0,
        )
        add_supplemental_attribute!(sys, br, o)
        push!(outages, o)
    end

    r_auto = solve_contingency_dc_power_flow(sys, outages; time_steps = 1, strategy = :auto)
    r_ref = solve_contingency_dc_power_flow(sys, outages; time_steps = 1, strategy = :refactor)

    for k in 2:size(r_auto.bus_angles, 3)
        @test isapprox(
            r_auto.bus_angles[:, :, k],
            r_ref.bus_angles[:, :, k];
            atol = 1e-10,
        )
        @test isapprox(
            r_auto.arc_active_power_flow_from_to[:, :, k],
            r_ref.arc_active_power_flow_from_to[:, :, k];
            atol = 1e-10,
        )
    end
end

@testset "Contingency DC PF — Woodbury N-2 matches rebuild" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    line1 = get_component(Line, sys, "1")
    line2 = get_component(Line, sys, "2")
    outage = GeometricDistributionForcedOutage(;
        mean_time_to_recovery = 0.0,
        outage_transition_probability = 0.0,
    )
    add_supplemental_attribute!(sys, line1, outage)
    add_supplemental_attribute!(sys, line2, outage)

    result =
        solve_contingency_dc_power_flow(sys, [outage]; time_steps = 1, strategy = :woodbury)
    @test all(result.converged)

    set_available!(line1, false)
    set_available!(line2, false)
    ref_data = PowerFlowData(DCPowerFlow(; time_steps = 1), sys)
    solve_power_flow!(ref_data)
    set_available!(line1, true)
    set_available!(line2, true)

    @test isapprox(result.bus_angles[:, :, 2], ref_data.bus_angles; atol = 1e-8)
end

@testset "Contingency DC PF — :auto_calibrate runs and matches :auto results" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    branches = collect(
        get_components(
            x -> !(typeof(x) <: Union{PhaseShiftingTransformer, DiscreteControlledACBranch}),
            ACTransmission,
            sys,
        ),
    )
    outages = PSY.Outage[]
    for br in branches
        o = GeometricDistributionForcedOutage(;
            mean_time_to_recovery = 0.0,
            outage_transition_probability = 0.0,
        )
        add_supplemental_attribute!(sys, br, o)
        push!(outages, o)
    end
    r_cal = solve_contingency_dc_power_flow(
        sys,
        outages;
        time_steps = 1,
        strategy = :auto_calibrate,
    )
    r_auto = solve_contingency_dc_power_flow(
        sys,
        outages;
        time_steps = 1,
        strategy = :auto,
    )
    @test all(r_cal.converged)
    @test isapprox(r_cal.bus_angles, r_auto.bus_angles; atol = 1e-10)
    @test isapprox(
        r_cal.arc_active_power_flow_from_to,
        r_auto.arc_active_power_flow_from_to;
        atol = 1e-10,
    )
end

@testset "Contingency DC PF — islanded and n_islands fields" begin
    # Non-islanding N-1 on c_sys5 (fully meshed) must report islanded=false
    # and n_islands=1 for every contingency.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    branches = collect(
        get_components(
            x -> !(typeof(x) <: Union{PhaseShiftingTransformer, DiscreteControlledACBranch}),
            ACTransmission,
            sys,
        ),
    )
    outages = PSY.Outage[]
    for br in branches
        o = GeometricDistributionForcedOutage(;
            mean_time_to_recovery = 0.0,
            outage_transition_probability = 0.0,
        )
        add_supplemental_attribute!(sys, br, o)
        push!(outages, o)
    end
    result = solve_contingency_dc_power_flow(sys, outages; time_steps = 1)
    for k in 2:length(PowerFlows.get_contingency_labels(result))
        @test PowerFlows.get_islanded(result, k) == false
        @test PowerFlows.get_n_islands(result, k) == 1
    end
    # Base slot defaults
    @test PowerFlows.get_islanded(result, 1) == false
    @test PowerFlows.get_n_islands(result, 1) == 1
end

@testset "Contingency DC PF — non-outaged arc flows match rebuild truth" begin
    # For each N-1 branch outage, confirm the per-arc flows stored in the
    # contingency result (at all arcs EXCEPT the outaged one) match the flows
    # obtained by rebuilding the system with set_available!(branch, false)
    # and running a fresh DC PF. This closes the gap left by angle-only
    # equivalence + outaged-arc-is-zero spot checks.
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # Snapshot the FULL-system arc axis BEFORE attaching outages or removing
    # branches. Used to build a permutation aligning rebuild-truth flows with
    # the contingency result's larger arc axis.
    snapshot_data = PowerFlowData(DCPowerFlow(; time_steps = 1), sys)
    orig_arcs = collect(PNM.get_arc_axis(snapshot_data.aux_network_matrix))
    orig_arc_lookup = Dict(a => i for (i, a) in enumerate(orig_arcs))

    branches = collect(
        get_components(
            x -> !(typeof(x) <: Union{PhaseShiftingTransformer, DiscreteControlledACBranch}),
            ACTransmission,
            sys,
        ),
    )
    outages = PSY.Outage[]
    for br in branches
        o = GeometricDistributionForcedOutage(;
            mean_time_to_recovery = 0.0,
            outage_transition_probability = 0.0,
        )
        add_supplemental_attribute!(sys, br, o)
        push!(outages, o)
    end

    result = solve_contingency_dc_power_flow(sys, outages; time_steps = 1)

    for (k0, br) in enumerate(branches)
        set_available!(br, false)
        ref = PowerFlowData(DCPowerFlow(; time_steps = 1), sys)
        solve_power_flow!(ref)
        ref_arcs = collect(PNM.get_arc_axis(ref.aux_network_matrix))
        set_available!(br, true)

        # Map each ref-system arc back to its position in the original axis.
        # The removed branch's arc is absent from ref_arcs by construction.
        perm = [orig_arc_lookup[a] for a in ref_arcs]

        flow_ctg = result.arc_active_power_flow_from_to[perm, 1, k0 + 1]
        flow_truth = ref.arc_active_power_flow_from_to[:, 1]
        @test isapprox(flow_ctg, flow_truth; atol = 1e-8)
    end
end

@testset "Contingency DC PF — parallel matches serial (determinism)" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    branches = collect(
        get_components(
            x -> !(typeof(x) <: Union{PhaseShiftingTransformer, DiscreteControlledACBranch}),
            ACTransmission,
            sys,
        ),
    )
    outages = PSY.Outage[]
    for br in branches
        o = GeometricDistributionForcedOutage(;
            mean_time_to_recovery = 0.0,
            outage_transition_probability = 0.0,
        )
        add_supplemental_attribute!(sys, br, o)
        push!(outages, o)
    end

    r_serial =
        solve_contingency_dc_power_flow(sys, outages; time_steps = 1, parallel = false)
    r_parallel =
        solve_contingency_dc_power_flow(sys, outages; time_steps = 1, parallel = true)

    @test r_serial.converged == r_parallel.converged
    @test isapprox(r_parallel.bus_angles, r_serial.bus_angles; atol = 1e-12)
    @test isapprox(
        r_parallel.arc_active_power_flow_from_to,
        r_serial.arc_active_power_flow_from_to;
        atol = 1e-12,
    )
end

@testset "Contingency DC PF — multi-timestep (angles match per-timestep rebuild)" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    line1 = get_component(Line, sys, "1")
    outage = GeometricDistributionForcedOutage(;
        mean_time_to_recovery = 0.0,
        outage_transition_probability = 0.0,
    )
    add_supplemental_attribute!(sys, line1, outage)

    pf = DCPowerFlow(; time_steps = 3)
    base_data = PowerFlowData(pf, sys)
    # Vary injections per timestep so solutions differ
    base_data.bus_active_power_injections[:, 2] .=
        base_data.bus_active_power_injections[:, 1] .* 1.05
    base_data.bus_active_power_injections[:, 3] .=
        base_data.bus_active_power_injections[:, 1] .* 0.95
    solve_power_flow!(base_data)

    n_buses = size(base_data.bus_angles, 1)
    n_arcs = size(base_data.arc_active_power_flow_from_to, 1)
    result = TimeContingencyPowerFlowData(
        n_buses,
        n_arcs,
        3,
        String["base", "out_line_1"];
        make_arc_active_power_losses = true,
    )
    solve_contingency_dc_power_flow!(result, base_data, sys, [outage])

    @test all(result.converged[:, 2])
    @test isapprox(result.bus_angles[:, :, 1], base_data.bus_angles; atol = 1e-10)
end
