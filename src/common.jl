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

function get_injections!(
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

function get_withdrawals!(
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

# sparse case (ABA and BA)

function my_mul_mt!(
    y::Vector{Float64},
    A::SparseMatrixCSC{Float64, Int64},
    x::Vector{Float64},
)
    copyto!(y, zeros(Float64, size(y)))
    for i in 1:size(A, 2)
        for j in A.colptr[i]:(A.colptr[i + 1] - 1)
            y[i] += A.nzval[j] * x[A.rowval[j]]
        end
    end
    return
end

# dense case (PTDF and ABA)

function my_mul_mt!(
    y::Vector{Float64},
    A::Matrix{Float64},
    x::Vector{Float64},
)
    y[:] .= transpose(A) * x
    return
end

# virtual case: all lines

function my_mul_mt!(
    y::Vector{Float64},
    A::PNM.VirtualPTDF,
    x::Vector{Float64},
)
    for i in eachindex(y)
        name_ = A.axes[1][i]
        y[i] = LinearAlgebra.dot(A[name_, :], x)
    end
    return
end

# virtual case: selected lines (not used yet)

function my_mul_mt!(
    y::Vector{Float64},
    A::PNM.VirtualPTDF,
    x::Vector{Float64},
    lines::Vector{String},
)
    for name_ in lines
        y[A.lookup[1][name_]] = LinearAlgebra.dot(A[name_, :], x)
    end
    return
end

# virtual case: single line (not used yet)

function my_mul_mt!(
    y::Vector{Float64},
    A::PNM.VirtualPTDF,
    x::Vector{Float64},
    line::String,
)
    y[A.lookup[1][line]] = LinearAlgebra.dot(A[line, :], x)
    return
end
