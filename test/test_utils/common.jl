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
