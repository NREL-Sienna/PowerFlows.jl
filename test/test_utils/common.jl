const SYSTEM_REIMPORT_COMPARISON_TOLERANCE = 1e-10
const POWERFLOW_COMPARISON_TOLERANCE = 1e-9

powerflow_match_fn(
    a::T,
    b::T,
) where {T <: Union{AbstractFloat, AbstractArray{<:AbstractFloat}}} =
    isapprox(a, b; atol = POWERFLOW_COMPARISON_TOLERANCE) || IS.isequivalent(a, b)
powerflow_match_fn(a, b) = IS.isequivalent(a, b)

# TODO temporary hack, see https://github.com/NREL-Sienna/PowerFlows.jl/issues/39
function PowerSystems.get_reactive_power_limits(gen::RenewableNonDispatch)
    gen_pf = get_power_factor(gen)
    gen_q = get_max_active_power(gen) * sqrt((1 / gen_pf^2) - 1)
    return (min = 0.0, max = gen_q)
end

# TODO more hacks
PowerSystems.get_r(::TwoTerminalHVDCLine) = 0.001
PowerSystems.get_x(::TwoTerminalHVDCLine) = 0.0
PowerSystems.get_b(::TwoTerminalHVDCLine) = (from = 0.0, to = 0.0)

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
    set_magnitude!(pq_bus, 0.54783)
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
    data.bus_magnitude[data.bus_lookup[117]] = 0.54783
    data.bus_angles[data.bus_lookup[117]] = 0.14956
end
