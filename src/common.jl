function get_total_p(l::PSY.PowerLoad)
    return PSY.get_active_power(l)
end

function get_total_q(l::PSY.PowerLoad)
    return PSY.get_reactive_power(l)
end

function get_total_p(l::PSY.StandardLoad)
    return PSY.get_constant_active_power(l) +
           PSY.get_current_active_power(l) +
           PSY.get_impedance_active_power(l)
end

function get_total_q(l::PSY.StandardLoad)
    return PSY.get_constant_reactive_power(l) +
           PSY.get_current_reactive_power(l) +
           PSY.get_impedance_reactive_power(l)
end

function get_total_p(l::PSY.ExponentialLoad)
    return PSY.get_active_power(l)
end

function get_total_q(l::PSY.ExponentialLoad)
    return PSY.get_reactive_power(l)
end

function _get_injections!(
    bus_activepower_injection::Vector{Float64},
    bus_reactivepower_injection::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System,
)
    sources = PSY.get_components(d -> !isa(d, PSY.ElectricLoad), PSY.StaticInjection, sys)
    for source in sources
        !PSY.get_available(source) && continue
        bus = PSY.get_bus(source)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        bus_activepower_injection[bus_ix] += PSY.get_active_power(source)
        bus_reactivepower_injection[bus_ix] += PSY.get_reactive_power(source)
    end
    return
end

function _get_withdrawals!(
    bus_activepower_withdrawals::Vector{Float64},
    bus_reactivepower_withdrawals::Vector{Float64},
    bus_lookup::Dict{Int, Int},
    sys::PSY.System,
)
    loads = PSY.get_components(x -> !isa(x, PSY.FixedAdmittance), PSY.ElectricLoad, sys)
    for l in loads
        !PSY.get_available(l) && continue
        bus = PSY.get_bus(l)
        bus_ix = bus_lookup[PSY.get_number(bus)]
        bus_activepower_withdrawals[bus_ix] += PSY.get_active_power(l)
        bus_reactivepower_withdrawals[bus_ix] += PSY.get_reactive_power(l)
    end
    return
end

##############################################################################
# Matrix Methods #############################################################

function my_mul_mt(
    A::PNM.VirtualPTDF,
    x::Vector{Float64},
)
    y = zeros(length(A.axes[1]))
    for i in 1:length(A.axes[1])
        name_ = A.axes[1][i]
        y[i] = LinearAlgebra.dot(A[name_, :], x)
    end
    return y
end