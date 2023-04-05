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
    @benchmark solve_powerflow!(DCPowerFlow(), pfd_1, pfd_1.power_network_matrix, pfd_1.aux_network_matrix, power_injections_1)
    @benchmark solve_powerflow!(DCPowerFlow(), pfd_1, pfd_1.power_network_matrix, pfd_1.aux_network_matrix, power_injections_1; parallel=true)
    solve_powerflow!(DCPowerFlow(), pfd_1, pfd_1.power_network_matrix, pfd_1.aux_network_matrix, power_injections_1)

    # with ABA
    @benchmark solve_powerflow!(DCPowerFlow(), pfd_2, pfd_2.power_network_matrix, pfd_2.aux_network_matrix, power_injections_2)
    solve_powerflow!(DCPowerFlow(), pfd_2, pfd_2.power_network_matrix, pfd_2.aux_network_matrix, power_injections_2)

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

"""
with no parallelization

* RTS *
- with ptdf
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  13.250 μs …  51.500 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     13.875 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   14.031 μs ± 946.094 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

     ▃▅█▇█▆▆▆▄▅▄▅▅▅▃▃▁▂▁                                        
  ▃▅▇████████████████████▇▇▇▆▅▆▇▇▄▄▃▃▄▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁ ▄
  13.2 μs         Histogram: frequency by time         15.8 μs <

 Memory estimate: 7.30 KiB, allocs estimate: 37.

- with ABA
BenchmarkTools.Trial: 10000 samples with 9 evaluations.
 Range (min … max):  2.199 μs … 735.606 μs  ┊ GC (min … max): 0.00% … 99.44%
 Time  (median):     2.370 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.608 μs ±  12.652 μs  ┊ GC (mean ± σ):  8.39% ±  1.72%

            ▂▃▅▆█▇▁███▇▇▅▅▄ ▂▁                                 
  ▂▁▂▂▃▃▄▅▅▇███████████████▆███▇▇▆▅▄▅▅▄▄▄▄▃▃▃▃▃▂▂▃▂▂▂▂▂▂▂▂▂▂▂ ▅
  2.2 μs          Histogram: frequency by time        2.71 μs <

 Memory estimate: 3.88 KiB, allocs estimate: 26.

 * ACTIVSg2000 *
 - with ptdf
 BenchmarkTools.Trial: 865 samples with 1 evaluation.
 Range (min … max):  5.680 ms …  6.343 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     5.759 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   5.768 ms ± 57.926 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

       ▁▃▅▆▅█▂▆▂▁▃▁▂▃▃▃▁                                      
  ▃▂▃▄▇█████████████████▆▇▆▄▄▃▂▂▃▃▃▁▂▂▁▂▁▁▁▁▁▁▂▁▂▂▁▁▂▁▁▁▁▁▁▂ ▄
  5.68 ms        Histogram: frequency by time        6.02 ms <

 Memory estimate: 133.51 KiB, allocs estimate: 44.

 - with ABA
 BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  31.792 μs …   8.317 ms  ┊ GC (min … max): 0.00% … 99.28%
 Time  (median):     33.458 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   38.113 μs ± 134.677 μs  ┊ GC (mean ± σ):  6.09% ±  1.72%

   ▄▇██▇▇▅▄▃▂▂▂▂▁▁          ▁▂▄▅▄▅▄▄▃▃▂▁▁▁▁▁▁▁▁                ▂
  ▆██████████████████▇▆█▇▇▇██████████████████████▇▇▆▅▅▅▆▆▆▅▄▅▃ █
  31.8 μs       Histogram: log(frequency) by time        49 μs <

 Memory estimate: 43.06 KiB, allocs estimate: 30.

 - with @threads
 BenchmarkTools.Trial: 252 samples with 1 evaluation.
 Range (min … max):  15.311 ms … 27.607 ms  ┊ GC (min … max):  0.00% … 20.76%
 Time  (median):     18.372 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   19.843 ms ±  4.201 ms  ┊ GC (mean ± σ):  11.68% ± 11.51%

  ▄▇█▂                                                         
  ████▆▄▄▃▃▂▃▃▂▃▃▃▃▄▄▃▃▂▁▃▃▁▁▁▂▁▁▁▃▄▁▂▄▃▃▃▃▄▃▄▃▄▇▅▆▄▆▅▂▄▃▂▄▁▃ ▃
  15.3 ms         Histogram: frequency by time        27.1 ms <

 Memory estimate: 49.42 MiB, allocs estimate: 3253.

* case_ACTIVSg10k *
- with ptdf
BenchmarkTools.Trial: 45 samples with 1 evaluation.
 Range (min … max):  111.699 ms … 115.260 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     112.024 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   112.275 ms ± 791.306 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

     ▄█▁▁▁▃                                                      
  ▄▆▇██████▄▆▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▄▁▁▁▁▁▄▁▁▁▁▁▄ ▁
  112 ms           Histogram: frequency by time          115 ms <

 Memory estimate: 650.80 KiB, allocs estimate: 47.

- with ABA
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  182.708 μs …  10.707 ms  ┊ GC (min … max): 0.00% … 97.82%
 Time  (median):     207.645 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   215.664 μs ± 123.973 μs  ┊ GC (mean ± σ):  0.76% ±  1.37%

   ▇█▆▅▃▃▃▄▄▁▁▁                                                  
  ▅███████████████▇███▇▇▇▇▆▆▆▅▅▄▄▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁ ▄
  183 μs           Histogram: frequency by time          306 μs <

 Memory estimate: 179.73 KiB, allocs estimate: 31.
"""

"""
with paralelization (Basic.Threads)

"""