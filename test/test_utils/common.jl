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

"Take RTS_GMLC_DA_sys and make some changes to it that are fully captured in the PowerFlowData(ACPowerFlow(), ...)"
function modify_rts_system!(sys::System)
    # For REF bus, voltage and angle are fixed; update active and reactive
    ref_bus = get_component(ACBus, sys, "Arne")  # bus number 113
    @assert get_bustype(ref_bus) == ACBusTypes.REF
    # NOTE: we are not testing the correctness of _power_redistribution_ref here, it is used on both sides of the test
    PF._power_redistribution_ref(sys, 2.4375, 0.1875, ref_bus)

    # For PV bus, active and voltage are fixed; update reactive and angle

    # For PQ bus, active and reactive are fixed; update voltage and angle

end

"Make the same changes to the PowerFlowData that modify_rts_system! makes to the System"
function modify_rts_powerflow!(data::PowerFlowData)
    data.bus_activepower_injection[data.bus_lookup[113]] -= 1.0
    data.bus_reactivepower_injection[data.bus_lookup[113]] -= 1.0
end
