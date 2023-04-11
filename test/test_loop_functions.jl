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
using SuiteSparse

const IS = InfrastructureSystems
const PSB = PowerSystemCaseBuilder
const PSY = PowerSystems
const PNM = PowerNetworkMatrices
const PF = PowerFlows

# get system

sys = System("ACTIVSg2000.m")
# sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

# use function
PowerFlowData(DCPowerFlow(), sys)

##############################################################################
# function step by step ######################################################

# first get the matrices
aux_network_matrix = PNM.BA_Matrix(sys)
power_network_matrix = PNM.ABA_Matrix(sys; factorize=true)

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

power_network_matrix.lookup[2].vals

# ! this does not work since there is a mismatch between lookup and temp_bus_map
for (ix, bus_no) in bus_lookup
    bus_name = temp_bus_map[bus_no]
    bus = PSY.get_component(PSY.Bus, sys, bus_name)
    bus_type[ix] = PSY.get_bustype(bus)
    bus_angle[ix] = PSY.get_angle(bus)
end

#! this is a problem also for "get_injections" and "get_withdrawals" functions

# consider lookup from the power_network_matrix
bus_activepower_injection = zeros(n_buses)          # !
bus_reactivepower_injection = zeros(n_buses)        # !
PF.get_injections!(bus_activepower_injection, bus_reactivepower_injection, bus_lookup, sys)

bus_activepower_withdrawals = zeros(n_buses)        # !
bus_reactivepower_withdrawals = zeros(n_buses)      # !
PF.get_withdrawals!(
    bus_activepower_withdrawals,
    bus_reactivepower_withdrawals,
    bus_lookup,
    sys,
)

##############################################################################
# test matrix multiplication: Sparse #########################################

# get the power injections
power_injection = bus_activepower_injection - bus_activepower_withdrawals

# get the angles
bus_angle = power_network_matrix.K \ power_injection[setdiff(1:end, power_network_matrix.ref_bus_positions)]

# get the flows
branch_flow_values = aux_network_matrix.data * bus_angle

# consider using for loop: my_transpose_mul_single!
y = zeros(length(branch_flow_values))
A = aux_network_matrix.data
x = bus_angle
for i in eachindex(y) # for each branch
    tmp = 0.0
    @show i
    for j in nzrange(A, i) # non zero bus indices
        @show A.rowval[j]
        tmp += A.nzval[j] * x[A.rowval[j]]
    end
    y[i] = tmp
end

# try to understand why...

i = 1                   # first branch --> select the first row of A
j_range = nzrange(A, i) # row range of non zeros for COLUMN i

# in fact I can see that...

isapprox([i for i in A[:, i] if i != 0], [A.nzval[k] for k in j_range], atol=1e-6)

# ...thus I am getting the column of A and not the row!

# IMPORTANT: nzrange(A, i) tells you the row range of the non zeros, not the indices of the corresponding nzval

# --> in fact...
i = 10
j_range = nzrange(A, i)
@assert [A.nzval[i] for i in j_range] == [i for i in A[:, i] if i != 0]

# now consider my_mul_single!
y = zeros(length(branch_flow_values))
A = aux_network_matrix.data
x = bus_angle
for i in eachindex(y)
    tmp = 0.0
    for j in A.colptr[i]:(A.colptr[i + 1] - 1) # non zero bus indices
        tmp += A.nzval[j] * x[A.rowval[j]]
    end
    y[i] = tmp
end

# try to understand why...
i = 1
j_range2 = A.colptr[i]:(A.colptr[i + 1] - 1)
j_range == j_range2 # A.colptr[i]:(A.colptr[i + 1] - 1) == nzrange(A, i)

# --> can be checked
for i in eachindex(y)
    jj1 = nzrange(A, i)
    jj2 = A.colptr[i]:(A.colptr[i + 1] - 1)
    if jj1 != jj2
        error("ranges are different")
    end
end


# SOLUTION --> finally try transpose 

A_t = sparse(transpose(A))

for i in eachindex(y) # for each branch
    tmp = 0.0
    for j in nzrange(A_t, i) # non zero bus indices
        tmp += A_t.nzval[j] * x[A_t.rowval[j]]
    end
    y[i] = tmp
end

# check
@assert isapprox(branch_flow_values, y, atol = 1e-6)