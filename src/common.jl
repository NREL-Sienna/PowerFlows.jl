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

################################## Matrix Methods ##########################################
function my_mul_mt!(
    y::Vector{Float64},
    A::SparseMatrixCSC{Float64, Int64},
    x::Vector{Float64},
)
    for i in 1:size(A, 2)
        row_idx = A.rowval[A.colptr[i]:(A.colptr[i + 1] - 1)]
        y[i] = LinearAlgebra.dot(A.nzval[A.colptr[i]:(A.colptr[i + 1] - 1)], x[row_idx])
    end
    return
end

function my_mul_mt_1!(
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

function my_mul_mt_2!(
    y::Vector{Float64},
    A::SparseMatrixCSC{Float64, Int64},
    x::Vector{Float64},
)
    y[:] .= transpose(A) * x
    return
end

# dense case #################################################################

function my_mul_mt!(
    y::Vector{Float64},
    A::Matrix{Float64},
    x::Vector{Float64}
)
    for i in eachindex(y)
        y[i] = LinearAlgebra.dot(A[:, i], x)
    end
    return
end

function my_mul_mt_1!(
    y::Vector{Float64},
    A::Matrix{Float64},
    x::Vector{Float64}
)
    copyto!(y, zeros(Float64, size(y)))
    for i in eachindex(y)
        temp = 0.0
        for j in eachindex(x)
            temp += A[j, i] * x[j]
        end
        y[i] = temp
    end
    return
end

function my_mul_mt_2!(
    y::Vector{Float64},
    A::Matrix{Float64},
    x::Vector{Float64},
)
    y[:] .= transpose(A) * x
    return
end