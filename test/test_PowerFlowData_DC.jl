using Revise
using Test
using Logging
using PowerSystems
using PowerSystemCaseBuilder
import PowerNetworkMatrices
using PowerFlows
using InfrastructureSystems
using LinearAlgebra
using CSV
using DataFrames
using BenchmarkTools

const IS = InfrastructureSystems
const PSB = PowerSystemCaseBuilder
const PSY = PowerSystems
const PNM = PowerNetworkMatrices

# @testset "test PowerFlowData structures" begin
# load test sytem

# sys = PSB.build_system(PSITestSystems, "test_RTS_GMLC_sys")
sys = PSY.System("ACTIVSg2000.m")
# sys = PSY.System("case_ACTIVSg10k.m")
# sys = System("case_ACTIVSg70k.m"; runchecks=false)

# get BA, ABA and PTDF matrix
BA = PNM.BA_Matrix(sys)
ABA = PNM.ABA_Matrix(sys; factorize = true)
ptdf = PNM.PTDF(sys)

# initialize PowerFlowData with PTDF matrix
# TODO: needs avaluation
pfd_1 = PowerFlowData(DCPowerFlow(), sys, ptdf, ABA)

# initialize PowerFlowData with ABA and BA matrices
# TODO: needs avaluation
pfd_2 = PowerFlowData(DCPowerFlow(), sys, ABA, BA)

# power_injections is the vector containing the solution of a UC and ED problem
# in this case we adopt the values in the pfd structure
power_injections_1 = pfd_1.bus_activepower_injection
power_injections_2 = pfd_2.bus_activepower_injection

# evaluate power flow

# with PTDF
@benchmark solve_powerflow!(
    DCPowerFlow(),
    pfd_1,
    pfd_1.power_network_matrix,
    pfd_1.aux_network_matrix,
    power_injections_1,
)
@benchmark solve_powerflow!(
    DCPowerFlow(),
    pfd_1,
    pfd_1.power_network_matrix,
    pfd_1.aux_network_matrix,
    power_injections_1;
    parallel = true,
)
solve_powerflow!(
    DCPowerFlow(),
    pfd_1,
    pfd_1.power_network_matrix,
    pfd_1.aux_network_matrix,
    power_injections_1,
)

# with ABA
@benchmark solve_powerflow!(
    DCPowerFlow(),
    pfd_2,
    pfd_2.power_network_matrix,
    pfd_2.aux_network_matrix,
    power_injections_2,
)
solve_powerflow!(
    DCPowerFlow(),
    pfd_2,
    pfd_2.power_network_matrix,
    pfd_2.aux_network_matrix,
    power_injections_2,
)

# check flows between the different methods
@test size(pfd_1.branch_flow_values) == size(pfd_2.branch_flow_values)

for br in 1:(pfd_1.n_branches)
    @test isapprox(
        pfd_1.branch_flow_values[br],
        pfd_2.branch_flow_values[br],
        atol = 1e-8,
    )
end

for b in 1:(pfd_1.n_buses - length(pfd_1.ref_buses))
    @test isapprox(pfd_1.bus_angle[b], pfd_2.bus_angle[b], atol = 1e-8)
end
# end
