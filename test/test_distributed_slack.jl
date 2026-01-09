
@testset "MULTI-PERIOD power flows evaluation with DS" begin
    for mode in (
            :nothing,
            :equal,
            :dict,
            :array_1,
            :array_24,
        ), ACSolver in (NewtonRaphsonACPowerFlow, TrustRegionACPowerFlow)
        @testset "AC Solver: $(ACSolver) mode: $(mode)" begin
            # get system
            sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
            generators = get_components(ThermalStandard, sys)

            create_gspf(generators, factor_values) =
                Dict(
                    (ThermalStandard, get_name(x)) => v for
                    (x, v) in zip(generators, factor_values)
                )
            Random.seed!(0)
            factor_values_1 = abs.(randn(Float64, length(generators)))
            Random.seed!(0)
            factor_values_24 = abs.(randn(Float64, length(generators), 24))
            generator_slack_participation_factors =
                if mode == :nothing
                    nothing
                elseif mode == :equal
                    Dict((ThermalStandard, get_name(x)) => 1.0 for x in generators)
                elseif mode == :dict
                    create_gspf(generators, factor_values_1)
                elseif mode == :array_1
                    [create_gspf(generators, factor_values_1)]
                elseif mode == :array_24
                    [create_gspf(generators, factor_values_24[:, t]) for t in 1:24]
                end

            # create structure for multi-period case
            time_steps = 24
            pf = ACPowerFlow(;
                generator_slack_participation_factors = generator_slack_participation_factors,
                time_steps = time_steps,
            )
            data = PowerFlowData(pf, sys)

            # allocate timeseries data from csv
            prepare_ts_data!(data, time_steps)

            init_p_injections = copy(data.bus_activepower_injection)
            solve_powerflow!(data)

            # check results
            subnetworks = PowerFlows._find_subnetworks_for_reference_buses(
                data.power_network_matrix.data,
                data.bus_type[:, 1],
            )
            for time_step in 1:time_steps
                _check_distributed_slack_consistency(
                    subnetworks,
                    data.bus_activepower_injection[:, time_step],
                    collect(data.bus_slack_participation_factors[:, time_step]),
                    init_p_injections[:, time_step],
                )
            end

            if mode == :array_24
                pf_short = ACPowerFlow(;
                    generator_slack_participation_factors = generator_slack_participation_factors[1:5],
                    time_steps = time_steps,
                )
                @test_throws ArgumentError(
                    "slack_participation_factors must have at least the same length as time_steps",
                ) PowerFlowData(pf_short, sys)
            end
        end
    end
end

function _get_spf_dict(sys::PSY.System,
    bus_numbers::Vector{Int},
    bus_slack_participation_factors::Vector{Float64},
)
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

@testset "AC PF with distributed slack" begin
    for (grid_lib, grid_name) in [
            (PSB.PSITestSystems, "c_sys14"),
            (PSB.MatpowerTestSystems, "matpower_case30_sys"),
        ], ACSolver in AC_SOLVERS_TO_TEST
        @testset "$(ACSolver) on $(grid_name)" begin
            sys = PSB.build_system(grid_lib, grid_name)
            # add a duplicate generator to a PV bus to make sure the method works for such set-ups
            g1 = first(
                get_components(
                    x -> get_bustype(get_bus(x)) == ACBusTypes.PV,
                    ThermalStandard,
                    sys,
                ),
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
                isempty(
                    get_components(x -> get_number(get_bus(x)) == bn, ThermalStandard, sys),
                ) &&
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

            pf = ACPowerFlow(; correct_bustypes = true)
            data = PowerFlowData(pf, sys)
            original_bus_power, original_gen_power =
                _system_generation_power(sys, bus_numbers)
            data_original_bus_power = copy(data.bus_activepower_injection[:, 1])
            res1 = solve_powerflow(pf, sys)

            bus_slack_participation_factors = zeros(Float64, length(bus_numbers))
            bus_slack_participation_factors[ref_n] .= 1.0

            pf2 = ACPowerFlow(;
                generator_slack_participation_factors = _get_spf_dict(
                    sys,
                    bus_numbers,
                    bus_slack_participation_factors,
                ),
                correct_bustypes = true,
            )
            res2 = solve_powerflow(pf2, sys)

            # basic test: if we pass the same slack participation factors as the default ones, the results
            # should be the same
            @test isapprox(
                res1["bus_results"].Vm,
                res2["bus_results"].Vm,
                atol = 1e-6,
                rtol = 0,
            )
            @test isapprox(
                res1["bus_results"].θ,
                res2["bus_results"].θ,
                atol = 1e-6,
                rtol = 0,
            )

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
                    sys,
                    bus_numbers,
                    bus_slack_participation_factors,
                ),
                correct_bustypes = true,
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
                    sys,
                    bus_numbers,
                    bus_slack_participation_factors,
                ),
                correct_bustypes = true,
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
                    sys,
                    bus_numbers,
                    bus_slack_participation_factors,
                ),
                correct_bustypes = true,
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
                        sys,
                        bus_numbers,
                        bus_slack_participation_factors,
                    ),
                ],
                correct_bustypes = true,
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
    end
end
