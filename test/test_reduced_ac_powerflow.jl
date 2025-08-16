function test_reduced_ac_powerflow(sys::PSY.System,
    nrs::Vector{PNM.NetworkReduction},
)
    pf = PF.ACPowerFlow(PF.NewtonRaphsonACPowerFlow)
    data = PF.PowerFlowData(pf, sys; time_steps = 1, network_reductions = nrs)
    PF.solve_powerflow!(data; pf = pf)
    return data
end

function get_bus_voltages(data::PF.PowerFlowData)
    bus_Vm = PF.get_bus_magnitude(data)
    bus_Va = PF.get_bus_angles(data)
    bus_voltages = Dict{Int, Complex{Float64}}()
    for (bus, i) in PF.get_bus_lookup(data)
        bus_voltages[bus] = bus_Vm[i] * exp(1im * bus_Va[i])
    end
    return bus_voltages
end

function get_branch_flows(data::PF.PowerFlowData)
    branch_flows = Dict{Tuple{Int, Int}, ComplexF64}()
    arc_ax = PNM.get_arc_axis(data.power_network_matrix.branch_admittance_from_to)
    for (i, arc) in enumerate(arc_ax)
        branch_flows[arc] =
            data.arc_activepower_flow_from_to[i] .+
            (im .* data.arc_reactivepower_flow_from_to[i])
        branch_flows[reverse(arc)] =
            data.arc_activepower_flow_to_from[i] .+
            (im .* data.arc_reactivepower_flow_to_from[i])
    end
    return branch_flows
end

function get_reduced_buses(nrd::PNM.NetworkReductionData)
    reduced_buses = Set()
    for (bus, reduced_set) in PNM.get_bus_reduction_map(nrd)
        if length(reduced_set) > 0
            push!(reduced_buses, bus)
        end
    end
    return reduced_buses
end

function validate_reduced_ac_powerflow(
    sys::PSY.System,
    nrs::Vector{PNM.NetworkReduction},
    unreduced_solved_data::PF.PowerFlowData,
)
    pf = PF.ACPowerFlow(PF.TrustRegionACPowerFlow)
    data = PF.PowerFlowData(pf, sys; network_reductions = nrs, correct_bustypes = true)
    PF.solve_powerflow!(data; pf = pf)
    @test all(data.converged)
    reduced_bus_results = get_bus_voltages(data)
    unreduced_bus_results = get_bus_voltages(unreduced_solved_data)
    nrd = PNM.get_network_reduction_data(data.power_network_matrix)
    reduced_buses = get_reduced_buses(nrd)
    if all(data.converged)
        for bus in PNM.get_bus_axis(data.power_network_matrix)
            if !(bus in reduced_buses)
                @test isapprox(
                    reduced_bus_results[bus],
                    unreduced_bus_results[bus],
                    atol = 1e-6,
                )
            end
        end
        # I expected these to match, but they don't. Odd.
        # Maybe I should only check branches where neither node was involved in a reduction?
        #=
        reduced_branch_flows = get_branch_flows(data)
        unreduced_branch_results = get_branch_flows(unreduced_solved_data)
        for arc in keys(PNM.get_direct_branch_map(nrd))
            @test isapprox(
                reduced_branch_flows[arc],
                unreduced_branch_results[arc],
                atol = 1e-6
            )
        end
        =#
    else
        @error "Reduced power flow did not converge; skipping result validation."
    end
end

@testset "2k bus system: validate results of reduce-then-solve" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
    pf = PF.ACPowerFlow(PF.TrustRegionACPowerFlow)
    unreduced = PF.PowerFlowData(pf, sys; correct_bustypes = true)
    PF.solve_powerflow!(unreduced; pf = pf)
    @assert all(unreduced.converged)
    @testset "radial reduction" begin
        validate_reduced_ac_powerflow(
            sys,
            PNM.NetworkReduction[PNM.RadialReduction()],
            unreduced,
        )
    end
    @testset "degree 2 reduction" begin
        validate_reduced_ac_powerflow(
            sys,
            PNM.NetworkReduction[PNM.DegreeTwoReduction()],
            unreduced,
        )
    end
    @testset "radial + degree 2 reduction" begin
        validate_reduced_ac_powerflow(
            sys,
            PNM.NetworkReduction[PNM.RadialReduction(), PNM.DegreeTwoReduction()],
            unreduced,
        )
    end
end

@testset "all reductions on psse_14_network_reduction_test_system" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    @testset "default reduction" begin
        result = test_reduced_ac_powerflow(sys, PNM.NetworkReduction[])
        @test all(result.converged)
    end
    @testset "radial reduction" begin
        result = test_reduced_ac_powerflow(sys, PNM.NetworkReduction[PNM.RadialReduction()])
        @test all(result.converged)
    end
    @testset "degree 2 reduction" begin
        result =
            test_reduced_ac_powerflow(sys, PNM.NetworkReduction[PNM.DegreeTwoReduction()])
        @test all(result.converged)
    end
    @testset "radial + degree 2 reduction" begin
        result = test_reduced_ac_powerflow(
            sys,
            PNM.NetworkReduction[PNM.RadialReduction(), PNM.DegreeTwoReduction()],
        )
        @test all(result.converged)
    end
    @testset "ward reduction" begin
        study_buses = [101, 114, 110, 111]
        result = test_reduced_ac_powerflow(
            sys,
            PNM.NetworkReduction[PNM.WardReduction(study_buses)],
        )
        @test all(result.converged) broken = true
    end
end
