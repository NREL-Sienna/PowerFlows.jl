const SYSTEM_REIMPORT_COMPARISON_TOLERANCE = 1e-10
const POWERFLOW_COMPARISON_TOLERANCE = 1e-9

powerflow_match_fn(
    a::T,
    b::T,
) where {T <: Union{AbstractFloat, AbstractArray{<:AbstractFloat}}} =
    isapprox(a, b; atol = POWERFLOW_COMPARISON_TOLERANCE) || IS.isequivalent(a, b)
powerflow_match_fn(a, b) = IS.isequivalent(a, b)

# TODO another temporary hack
"Create a version of the RTS_GMLC system that plays nice with the current implementation of AC power flow"
function create_pf_friendly_rts_gmlc()
    sys = build_system(PSISystems, "RTS_GMLC_DA_sys")
    remove_component!(sys, only(get_components(PSY.TwoTerminalHVDC, sys)))  # HVDC power flow not implemented yet
    # Modify some things so reactive power redistribution succeeds
    for (component_type, component_name, new_limits) in [
        (RenewableDispatch, "113_PV_1", (min = -30.0, max = 30.0))
        (ThermalStandard, "115_STEAM_3", (min = -50.0, max = 100.0))
        (ThermalStandard, "207_CT_1", (min = -70.0, max = 70.0))
        (RenewableDispatch, "215_PV_1", (min = -40.0, max = 40.0))
        (ThermalStandard, "307_CT_1", (min = -70.0, max = 70.0))
        (ThermalStandard, "315_CT_8", (min = 0.0, max = 80.0))
    ]
        set_reactive_power_limits!(
            get_component(component_type, sys, component_name),
            new_limits,
        )
    end

    # Patch https://github.com/NREL-Sienna/PowerFlows.jl/issues/47
    sync_conds = filter(c -> occursin("SYNC_COND", get_name(c)),
        collect(get_components(StaticInjection, sys)))
    set_base_power!.(sync_conds, 100.0)
    return sys
end

"Take RTS_GMLC_DA_sys and make some changes to it that are fully captured in the PowerFlowData(ACPowerFlow(), ...)"
function modify_rts_system!(sys::System)
    # For REF bus, voltage and angle are fixed; update active and reactive
    ref_bus = get_bus(sys, 113)  # "Arne"
    @assert get_bustype(ref_bus) == ACBusTypes.REF
    # NOTE: we are not testing the correctness of _power_redistribution_ref here, it is used on both sides of the test
    PF._power_redistribution_ref(
        sys,
        2.4375,
        0.1875,
        ref_bus,
        PF.DEFAULT_MAX_REDISTRIBUTION_ITERATIONS,
    )

    # For PV bus, active and voltage are fixed; update reactive and angle
    pv_bus = get_bus(sys, 202)  # "Bacon"
    @assert get_bustype(pv_bus) == ACBusTypes.PV
    PF._reactive_power_redistribution_pv(
        sys,
        0.37267,
        pv_bus,
        PF.DEFAULT_MAX_REDISTRIBUTION_ITERATIONS,
    )
    set_angle!(pv_bus, -0.13778)

    # For PQ bus, active and reactive are fixed; update voltage and angle
    pq_bus = get_bus(sys, 117)  # "Aston"
    @assert get_bustype(pq_bus) == ACBusTypes.PQ
    set_magnitude!(pq_bus, 0.84783)
    set_angle!(pq_bus, 0.14956)
end

"Make the same changes to the PowerFlowData that modify_rts_system! makes to the System"
function modify_rts_powerflow!(data::PowerFlowData)
    # For REF bus, voltage and angle are fixed; update active and reactive
    data.bus_activepower_injection[data.bus_lookup[113]] = 2.4375
    data.bus_reactivepower_injection[data.bus_lookup[113]] = 0.1875

    # For PV bus, active and voltage are fixed; update reactive and angle
    data.bus_reactivepower_injection[data.bus_lookup[202]] = 0.37267
    data.bus_angles[data.bus_lookup[202]] = -0.13778

    # For PQ bus, active and reactive are fixed; update voltage and angle
    data.bus_magnitude[data.bus_lookup[117]] = 0.84783
    data.bus_angles[data.bus_lookup[117]] = 0.14956
end

function _system_generation_power(
    sys::System,
    bus_numbers::Vector{Int},
)
    bus_power = zeros(Float64, length(bus_numbers))
    generators = collect(get_components(Union{Generator, Source}, sys))
    gen_power = zeros(Float64, length(generators))
    with_units_base(sys, UnitSystem.NATURAL_UNITS) do
        bus_power .= [
            isempty(g) ? 0 : sum([get_active_power(gg) for gg in g]) for g in [
                get_components(
                    x -> get_number(get_bus(x)) == i,
                    Union{Generator, Source},
                    sys,
                )
                for i in bus_numbers
            ]
        ]
        gen_power .= get_active_power.(generators)
    end
    return bus_power, gen_power
end

function _reset_gen_power!(
    sys::System,
    original_gen_power::Vector{Float64},
)
    with_units_base(sys, UnitSystem.NATURAL_UNITS) do
        for (g, og) in
            zip(get_components(Union{Generator, Source}, sys), original_gen_power)
            set_active_power!(g, og)
        end
    end
end

function _check_distributed_slack_consistency(
    subnetworks::Dict{Int, Vector{Int}},
    result_bus_power::Vector{Float64},
    slack_participation_factors::Vector{Float64},
    original_bus_power::Vector{Float64},
)
    for (_, subnetwork_buses) in subnetworks
        subnetwork_factors = slack_participation_factors[subnetwork_buses]
        slack_provided =
            result_bus_power[subnetwork_buses] .- original_bus_power[subnetwork_buses]
        nnz = subnetwork_factors .!= 0.0

        @test all(isapprox.(slack_provided[.!nnz], 0.0, atol = 1e-6, rtol = 0))
        @test !any(isapprox.(slack_provided[nnz], 0.0, atol = 1e-6, rtol = 0))
        @test isapprox(
            slack_provided[nnz] ./ sum(slack_provided),
            subnetwork_factors[nnz] ./ sum(subnetwork_factors);
            atol = 1e-6,
            rtol = 0,
        )
    end
    return
end

function _check_ds_pf(
    pf::ACPowerFlow,
    sys::System,
    bus_slack_participation_factors::Vector{Float64},
    bus_numbers::Vector{Int},
    original_bus_power::Vector{Float64},
    original_gen_power::Vector{Float64},
    data_original_bus_power::Vector{Float64};
    kwargs...,
)
    res = solve_powerflow(pf, sys; kwargs...)

    data = PowerFlowData(pf, sys; kwargs...)
    subnetworks = PowerFlows._find_subnetworks_for_reference_buses(
        data.power_network_matrix.data,
        data.bus_type[:, 1],
    )

    _check_distributed_slack_consistency(
        subnetworks,
        res["bus_results"][:, :P_gen],
        bus_slack_participation_factors,
        original_bus_power,
    )

    solve_powerflow!(pf, sys; kwargs...)
    p_solve, _ = _system_generation_power(sys, bus_numbers)

    @test isapprox(p_solve, res["bus_results"][:, :P_gen]; atol = 1e-6, rtol = 0)

    _reset_gen_power!(sys, original_gen_power)
    # to make sure the reset function is working properly:
    p_bus_reset, p_gen_reset = _system_generation_power(sys, bus_numbers)
    @test original_bus_power == p_bus_reset
    @test original_gen_power == p_gen_reset

    @test data.bus_slack_participation_factors[:, 1] == bus_slack_participation_factors
    solve_powerflow!(data; pf = pf)
    # now check the slack power distribution logic
    _check_distributed_slack_consistency(
        subnetworks,
        data.bus_activepower_injection[:, 1],
        bus_slack_participation_factors,
        data_original_bus_power,
    )
    return
end

"""These functions are used to create simple components for the tests to have more compact code"""

"""
    _check_name(sys::System, name::String, component_type::DataType)
    Check if the name is unique in the system. If not, append a number to the name.
"""
function _check_name(sys::System, name::String, component_type::DataType)
    # Check if the name is unique
    check = true
    i = 1
    while check
        if has_component(sys, component_type, name)
            i += 1
            name = name * "_$i"
        else
            check = false
        end
    end
    return name
end

"""    
    _add_simple_bus!(sys::System, number::Int, bus_type::ACBusTypes, base_voltage::Number, voltage_magnitude::Float64=1.0, voltage_angle::Float64=0.0)
    Simplified function to create and add a bus to the system with the given parameters.
"""
function _add_simple_bus!(
    sys::System,
    number::Int,
    bus_type::ACBusTypes,
    base_voltage::Number,
    voltage_magnitude::Float64 = 1.0,
    voltage_angle::Float64 = 0.0,
)
    bus = ACBus(;
        number = number,
        name = _check_name(sys, "bus_$number", ACBus),
        available = true,
        bustype = bus_type,
        angle = voltage_angle,
        magnitude = voltage_magnitude,
        voltage_limits = (0.0, 2.0),
        base_voltage = Float64(base_voltage),
    )
    add_component!(sys, bus)
    return bus
end

"""    
    _add_simple_load!(sys::System, bus::ACBus, active_power::Number, reactive_power::Number)
    Simplified function to create and add a load to the system with the given parameters.
"""
function _add_simple_load!(
    sys::System,
    bus::ACBus,
    active_power::Number,
    reactive_power::Number,
)
    load = PowerLoad(;
        name = _check_name(sys, "load_$(get_number(bus))", PowerLoad),
        available = true,
        bus = bus,
        active_power = Float64(active_power), # Per-unitized by device base_power
        reactive_power = Float64(reactive_power), # Per-unitized by device base_power
        base_power = 1.0, # MVA
        max_active_power = 100.0, # 10 MW per-unitized by device base_power
        max_reactive_power = 100.0,
    )

    add_component!(sys, load)
    return load
end

"""    
    _add_simple_source!(sys::System, bus::ACBus, active_power::Number=0.0, reactive_power::Number=0.0)
    Simplified function to create and add a source to the system with the given parameters.
"""
function _add_simple_source!(
    sys::System,
    bus::ACBus,
    active_power::Number = 0.0,
    reactive_power::Number = 0.0,
)
    source = Source(;
        name = _check_name(sys, "source_$(get_number(bus))", Source),
        available = true,
        bus = bus,
        active_power = Float64(active_power),
        reactive_power = Float64(reactive_power),
        R_th = 1e-5,
        X_th = 1e-5,
    )
    add_component!(sys, source)
    return source
end

"""    
    _add_simple_thermal_standard!(sys::System, bus::ACBus, active_power::Number=0.0, reactive_power::Number=0.0)
    Simplified function to create and add a thermal standard generator to the system with the given parameters.
"""
function _add_simple_thermal_standard!(
    sys::System,
    bus::ACBus,
    active_power::Number,
    reactive_power::Number,
)
    gen = ThermalStandard(;
        name = _check_name(sys, "thermal_standard_$(get_number(bus))", ThermalStandard),
        available = true,
        status = true,
        bus = bus,
        active_power = Float64(active_power),
        reactive_power = Float64(reactive_power),
        rating = 1.0,
        active_power_limits = (min = 0, max = 1),
        reactive_power_limits = (min = -1, max = 1),
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
    add_component!(sys, gen)
    return gen
end

"""    
    _add_simple_line!(sys::System, bus1::ACBus, bus2::ACBus, r::Float64=1e-3, x::Float64=1e-3, b::Float64=0.0)
    Simplified function to create and add a line to the system with the given parameters.
"""
function _add_simple_line!(
    sys::System,
    bus1::ACBus,
    bus2::ACBus,
    r::Float64 = 1e-3,
    x::Float64 = 1e-3,
    b::Float64 = 0.0,
)
    line = Line(;
        name = _check_name(sys, "line_$(get_number(bus1))_$(get_number(bus2))", Line),
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = Arc(; from = bus1, to = bus2),
        r = r,
        x = x,
        b = (from = b / 2, to = b / 2),
        rating = 1.0,
        angle_limits = (min = -pi / 2, max = pi / 2),
    )
    add_component!(sys, line)
    return line
end
