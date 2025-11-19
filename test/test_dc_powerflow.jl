function flows_from_dataframe(flow_results_df::DataFrame,
    arc_lookup::Dict{Tuple{Int, Int}, Int},
    direction::Symbol = :P_from_to,
)
    flows = fill(NaN, length(arc_lookup))
    for row in eachrow(flow_results_df)
        flows[arc_lookup[(row.bus_from, row.bus_to)]] = row[direction]
    end
    @assert !any(isnan.(flows))
    return flows
end

@testset "SINGLE PERIOD power flows evaluation: ABA, PTDF, VirtualPTDF" begin
    # get system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

    # get indices
    buses = collect(PSY.get_components(PSY.ACBus, sys))

    # get sorted indices for branches
    branches = collect(PSY.get_available_components(PSY.ACBranch, sys))
    from_bus = [PSY.get_number(PSY.get_arc(x).from) for x in branches]

    # get reference values: flows and angles.
    # See issue 210: would be better to compare against external program.
    data = PowerFlowData(DCPowerFlow(), sys; correct_bustypes = true)
    power_injection =
        deepcopy(data.bus_activepower_injection - data.bus_activepower_withdrawals)
    matrix_data = deepcopy(data.power_network_matrix.K)       # LU factorization of ABA
    aux_network_matrix = deepcopy(data.aux_network_matrix)    # BA matrix

    valid_ix = setdiff(
        1:length(power_injection),
        PNM.get_ref_bus_position(data.aux_network_matrix),
    )
    ref_bus_angles = deepcopy(data.bus_angles)
    ref_flow_values = deepcopy(data.arc_activepower_flow_from_to)

    ref_bus_angles[valid_ix] = matrix_data \ power_injection[valid_ix]
    ref_flow_values = transpose(aux_network_matrix.data) * ref_bus_angles

    basepower = PSY.get_base_power(sys)
    arc_lookup = PF.get_arc_lookup(data)
    # CASE 1: ABA and BA matrices
    solved_data_ABA = solve_powerflow(DCPowerFlow(), sys; correct_bustypes = true)
    ABA_branch_flows = solved_data_ABA["1"]["flow_results"]
    @test isapprox(
        1 / basepower .* flows_from_dataframe(ABA_branch_flows, arc_lookup, :P_from_to),
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        1 / basepower .* flows_from_dataframe(ABA_branch_flows, arc_lookup, :P_to_from),
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_ABA["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)

    # CASE 2: PTDF and ABA MATRICES
    solved_data_PTDF = solve_powerflow(PTDFDCPowerFlow(), sys; correct_bustypes = true)
    PTDF_branch_flows = solved_data_PTDF["1"]["flow_results"]
    @test isapprox(
        1 / basepower .* flows_from_dataframe(PTDF_branch_flows, arc_lookup, :P_from_to),
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        1 / basepower .* flows_from_dataframe(PTDF_branch_flows, arc_lookup, :P_to_from),
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_PTDF["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)

    # CASE 3: VirtualPTDF and ABA MATRICES
    solved_data_vPTDF = solve_powerflow(vPTDFDCPowerFlow(), sys; correct_bustypes = true)
    vPTDF_branch_flows = solved_data_vPTDF["1"]["flow_results"]
    @test isapprox(
        1 / basepower .* flows_from_dataframe(vPTDF_branch_flows, arc_lookup, :P_from_to),
        ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(
        1 / basepower .* flows_from_dataframe(vPTDF_branch_flows, arc_lookup, :P_to_from),
        -ref_flow_values,
        atol = 1e-6,
    )
    @test isapprox(solved_data_vPTDF["1"]["bus_results"].θ, ref_bus_angles, atol = 1e-6)
end

@testset "DC power flow with an LCC" begin
    sys, _ = simple_lcc_system()
    lcc = first(get_components(PSY.TwoTerminalLCCLine, sys))
    with_units_base(sys, PSY.UnitSystem.SYSTEM_BASE) do
        set_active_power_flow!(lcc, 0.3)
    end
    for T in (DCPowerFlow, PTDFDCPowerFlow, vPTDFDCPowerFlow)
        results = solve_powerflow(T(), sys; correct_bustypes = true)
        lcc_flow = results["1"]["lcc_results"][1, :P_from_to]
        with_units_base(sys, PSY.UnitSystem.NATURAL_UNITS) do
            @test lcc_flow == get_active_power_flow(lcc)
        end
    end
end

# TODO LCC DC test case with nonzero loss.

@testset "DC power flow: results independent of units" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    for T in (DCPowerFlow, PTDFDCPowerFlow, vPTDFDCPowerFlow)
        line_name, flow_natural =
            power_flow_with_units(sys, T, PSY.UnitSystem.NATURAL_UNITS)
        line_name2, flow_system = power_flow_with_units(sys, T, PSY.UnitSystem.SYSTEM_BASE)
        @test line_name == line_name2
        @test flow_natural == flow_system
    end
end
