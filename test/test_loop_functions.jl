using Test
using Logging
using PowerSystems
using PowerSystemCaseBuilder
using PowerNetworkMatrices
using PowerFlows
using InfrastructureSystems
using LinearAlgebra
using CSV
using DataFrames
using SparseArrays
using Plots

const IS = InfrastructureSystems
const PSB = PowerSystemCaseBuilder
const PSY = PowerSystems
const PNM = PowerNetworkMatrices

function my_transpose_mul_la!(y::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64})
    y[:] .= transpose(A)*x
    return
end

function my_transpose_mul_single!(y::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64})
    for i in eachindex(y) # for each branch
        tmp = 0.0
        for j in nzrange(A, i) # non zero bus indices
            tmp += A.nzval[j] * x[A.rowval[j]]
        end
        y[i] = tmp
    end
    return
end

function my_transpose_mul_mt!(y::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64})
    Threads.@threads for i in eachindex(y) # for each branch
        tmp = 0.0
        for j in nzrange(A, i) # non zero bus indices
            tmp += A.nzval[j] * x[A.rowval[j]]
        end
        y[i] = tmp
    end
    return
end

function my_mul_la!(y::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64})
    mul!(y, A, x)
    return
end

function my_mul_single!(y::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64})
    for i in eachindex(y) # for each bus
        tmp = 0.0
        for j in A.colptr[i]:(A.colptr[i+1]-1) # non zero bus indices
            tmp += A.nzval[j] * x[A.rowval[j]]
        end
        y[i] = tmp
    end
    return
end

function my_mul_mt!(y::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64})
    Threads.@threads for i in eachindex(y) # for each bus
        tmp = 0.0
        for j in A.colptr[i]:(A.colptr[i+1]-1) # non zero bus indices
            tmp += A.nzval[j] * x[A.rowval[j]]
        end
        y[i] = tmp
    end
    return
end

# get system

sys = System("ACTIVSg2000.m")

# use function
PowerFlowData(DCPowerFlow(), sys)

##############################################################################
# function step by step ######################################################

# first get the matrices
power_network_matrix = PNM.ABA_Matrix(sys)
aux_network_matrix = PNM.BA_Matrix(sys)

# check if the maps betwen the 2 matrices match

# this is the way values are evaluated in BA_Matrix
line_ax = [PSY.get_name(branch) for branch in branches]
bus_ax = [PSY.get_number(bus) for bus in setdiff(buses, aux_network_matrix.ref_bus_positions)]
lookup = (PNM.make_ax_ref(line_ax), PNM.make_ax_ref(bus_ax))
# ! they are different
@assert length(setdiff(aux_network_matrix.lookup[1].vals, lookup[1].vals)) == 0
@assert length(setdiff(aux_network_matrix.lookup[2].vals, lookup[2].vals)) == 0

# and this is how values are evaluated in ABA_Matrix
line_ax1 = [PSY.get_name(branch) for branch in branches]
bus_ax1 = [PSY.get_number(bus) for bus in buses]
axes1 = (line_ax, setdiff(bus_ax, power_network_matrix.ref_bus_positions))
lookup1 = (PNM.make_ax_ref(line_ax), PNM.make_ax_ref(bus_ax))
# ! they are different
@assert length(setdiff(power_network_matrix.lookup[1].vals, lookup[1].vals)) == 0
@assert length(setdiff(power_network_matrix.lookup[2].vals, lookup[2].vals)) == 0

# get number of buses and branches
n_buses = length(axes(power_network_matrix.data, 2)) + length(power_network_matrix.ref_bus_positions)
n_branches = length(axes(aux_network_matrix.data, 1))

bus_lookup = power_network_matrix.lookup[2]
branch_lookup = power_network_matrix.lookup[1]
bus_type = Vector{Any}(undef, n_buses) # error with PSY.BusTypes
bus_angle = Vector{Float64}(undef, n_buses)
temp_bus_map = Dict{Int, String}(
    PSY.get_number(b) => PSY.get_name(b) for b in PSY.get_components(PSY.Bus, sys)
)

# ! this does not work since there is a mismatch between lookup and dictionary
for (ix, bus_no) in bus_lookup
    bus_name = temp_bus_map[bus_no]
    bus = PSY.get_component(PSY.Bus, sys, bus_name)
    bus_type[ix] = PSY.get_bustype(bus)
    bus_angle[ix] = PSY.get_angle(bus)
end

#! this is a problem also for "get_injections" and "get_withdrawals" functions

##############################################################################
# test matrix multiplication: Sparse #########################################

# # get the power injections
# power_injection = # ! to be done

# # get the angles
# data.bus_angle[:] .= matrix_data \ power_injection

# # get the flows
# data.branch_flow_values[:] = aux_network_matrix.data * data.bus_angle