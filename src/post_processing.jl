
"""
Calculates the From - To complex power flow (Flow injected at the bus) of branch of type
TapTransformer
"""
function flow_val(b::TapTransformer)
    !get_available(b) && return 0.0
    Y_t = get_series_admittance(b)
    c = 1 / get_tap(b)
    arc = get_arc(b)
    V_from = arc.from.magnitude * (cos(arc.from.angle) + sin(arc.from.angle) * 1im)
    V_to = arc.to.magnitude * (cos(arc.to.angle) + sin(arc.to.angle) * 1im)
    I = (V_from * Y_t * c^2) - (V_to * Y_t * c)
    flow = V_from * conj(I)
    return flow
end

"""
Calculates the From - To complex power flow (Flow injected at the bus) of branch of type
Line
"""
function flow_val(b::ACBranch)
    !get_available(b) && return 0.0
    Y_t = get_series_admittance(b)
    arc = PowerSystems.get_arc(b)
    V_from = arc.from.magnitude * (cos(arc.from.angle) + sin(arc.from.angle) * 1im)
    V_to = arc.to.magnitude * (cos(arc.to.angle) + sin(arc.to.angle) * 1im)
    I = V_from * (Y_t + (1im * PowerSystems.get_b(b).from)) - V_to * Y_t
    flow = V_from * conj(I)
    return flow
end

"""
Calculates the From - To complex power flow (Flow injected at the bus) of branch of type
Transformer2W
"""
function flow_val(b::Transformer2W)
    !get_available(b) && return 0.0
    Y_t = get_series_admittance(b)
    arc = get_arc(b)
    V_from = arc.from.magnitude * (cos(arc.from.angle) + sin(arc.from.angle) * 1im)
    V_to = arc.to.magnitude * (cos(arc.to.angle) + sin(arc.to.angle) * 1im)
    I = V_from * (Y_t + (1im * get_primary_shunt(b))) - V_to * Y_t
    flow = V_from * conj(I)
    return flow
end

function flow_val(b::PhaseShiftingTransformer)
    error("Systems with PhaseShiftingTransformer not supported yet")
    return
end

"""
Calculates the From - To complex power flow using external data of voltages of branch of type
TapTransformer
"""
function flow_func(b::TapTransformer, V_from::Complex{Float64}, V_to::Complex{Float64})
    !get_available(b) && return (0.0, 0.0)
    Y_t = get_series_admittance(b)
    c = 1 / get_tap(b)
    I = (V_from * Y_t * c^2) - (V_to * Y_t * c)
    flow = V_from * conj(I)
    return real(flow), imag(flow)
end

"""
Calculates the From - To complex power flow using external data of voltages of branch of type
Line
"""
function flow_func(b::ACBranch, V_from::Complex{Float64}, V_to::Complex{Float64})
    !get_available(b) && return (0.0, 0.0)
    Y_t = get_series_admittance(b)
    I = V_from * (Y_t + (1im * get_b(b).from)) - V_to * Y_t
    flow = V_from * conj(I)
    return real(flow), imag(flow)
end

"""
Calculates the From - To complex power flow using external data of voltages of branch of type
Transformer2W
"""
function flow_func(b::Transformer2W, V_from::Complex{Float64}, V_to::Complex{Float64})
    !get_available(b) && return (0.0, 0.0)
    Y_t = get_series_admittance(b)
    I = V_from * (Y_t + (1im * get_primary_shunt(b))) - V_to * Y_t
    flow = V_from * conj(I)
    return real(flow), imag(flow)
end

function flow_func(
    b::PhaseShiftingTransformer,
    V_from::Complex{Float64},
    V_to::Complex{Float64},
)
    error("Systems with PhaseShiftingTransformer not supported yet")
    return
end

"""
Updates the flow on the branches
"""
function _update_branch_flow!(sys::System)
    for b in get_components(ACBranch, sys)
        S_flow = get_available(b) ? flow_val(b) : 0.0 + 0.0im
        set_active_power_flow!(b, real(S_flow))
        set_reactive_power_flow!(b, imag(S_flow))
    end
end

"""
Obtain total load on bus b
"""
function _get_load_data(sys::System, b::PSY.Bus)
    active_power = 0.0
    reactive_power = 0.0
    for l in get_components(ElectricLoad, sys, x -> !isa(x, FixedAdmittance))
        !get_available(l) && continue
        if (l.bus == b)
            active_power += get_active_power(l)
            reactive_power += get_reactive_power(l)
        end
    end
    return active_power, reactive_power
end

function _get_fixed_admittance_power(sys::System, b::PSY.Bus, result::AbstractVector, ix::Int)
    active_power = 0.0
    reactive_power = 0.0
    for l in get_components(FixedAdmittance, sys)
        !get_available(l) && continue
        if (l.bus == b)
            Vm_squared =
                b.bustype == PSY.BusTypes.PQ ? result[2 * ix - 1]^2 : get_magnitude(b)^2
            active_power += Vm_squared * real(get_Y(l))
            reactive_power += Vm_squared * imag(get_Y(l))
        end
    end
    return active_power, reactive_power
end
