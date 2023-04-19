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
using BenchmarkTools

const IS = InfrastructureSystems
const PSB = PowerSystemCaseBuilder
const PSY = PowerSystems
const PNM = PowerNetworkMatrices
const PF = PowerFlows

# get system

# sys = System("ACTIVSg2000.m")
# sys = System("case_ACTIVSg10k.m")
sys = System("case_ACTIVSg70k.m"; runchecks=false)
# sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)

##############################################################################
# CASE 1: ABA and BA matrices ################################################
data = PowerFlowData(DCPowerFlow(), sys)
data_1 = PowerFlowData(DCPowerFlow(), sys)
data_2 = PowerFlowData(DCPowerFlow(), sys)
power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
matrix_data = data.power_network_matrix.K       # LU factorization of ABA
aux_network_matrix = data.aux_network_matrix    # BA matrix
valid_ix = setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
data.bus_angle[valid_ix] = matrix_data \ power_injection[valid_ix]

# comparison
@benchmark PF.my_mul_mt!(data.branch_flow_values, aux_network_matrix.data, data.bus_angle)
@benchmark PF.my_mul_mt_1!(data_1.branch_flow_values, aux_network_matrix.data, data.bus_angle)
@benchmark PF.my_mul_mt_2!(data_2.branch_flow_values, aux_network_matrix.data, data.bus_angle)

PF.my_mul_mt!(data.branch_flow_values, aux_network_matrix.data, data.bus_angle)
PF.my_mul_mt_1!(data_1.branch_flow_values, aux_network_matrix.data, data.bus_angle)
PF.my_mul_mt_2!(data_2.branch_flow_values, aux_network_matrix.data, data.bus_angle)

isapprox(data.branch_flow_values, data_2.branch_flow_values, atol = 1e-5)
isapprox(data.branch_flow_values, data_1.branch_flow_values, atol = 1e-6)
isapprox(data_1.branch_flow_values, data_2.branch_flow_values, atol = 1e-6)
isapprox(data_1.branch_flow_values, data.branch_flow_values, atol = 1e-6)

sum(abs.(data.branch_flow_values - data_2.branch_flow_values))
sum(abs.(data.branch_flow_values - data_1.branch_flow_values))
sum(abs.(data_2.branch_flow_values - data_1.branch_flow_values))

"""
2k
BenchmarkTools.Trial: 9657 samples with 1 evaluation.
 Range (min … max):  457.083 μs …   8.791 ms  ┊ GC (min … max): 0.00% … 92.32%
 Time  (median):     471.500 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   516.253 μs ± 591.795 μs  ┊ GC (mean ± σ):  8.18% ±  6.71%

    ▂▅▇▅▂█▄  ▂▅▆▅▃                                               
  ▂▄███████████████▆▅▄▃▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  457 μs           Histogram: frequency by time          540 μs <

 Memory estimate: 752.05 KiB, allocs estimate: 9622.

 BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  11.750 μs …  6.811 ms  ┊ GC (min … max): 0.00% … 99.21%
 Time  (median):     12.500 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   14.217 μs ± 68.014 μs  ┊ GC (mean ± σ):  4.75% ±  0.99%

    ▅██▆▄                                                      
  ▁▄██████▅▄▃▃▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▂▃▃▃▃▂▃▂▃▃▂▂▂▂▁▁▁▁▁▁ ▂
  11.8 μs         Histogram: frequency by time          19 μs <

 Memory estimate: 25.75 KiB, allocs estimate: 6.

 BenchmarkTools.Trial: 10000 samples with 4 evaluations.
 Range (min … max):   7.083 μs …  1.733 ms  ┊ GC (min … max):  0.00% … 98.88%
 Time  (median):      9.396 μs              ┊ GC (median):     0.00%
 Time  (mean ± σ):   11.256 μs ± 44.899 μs  ┊ GC (mean ± σ):  10.49% ±  2.62%

  ▂▇█▇▅▃▁        ▁▂▁▁               ▄▆▇▇▆▅▅▃▂▁▁               ▃
  ████████████████████▇▄▆▅▁▅▅▅▅▁▃▃▄████████████████▇▆▆▇▆▅▁▅▇▇ █
  7.08 μs      Histogram: log(frequency) by time      15.6 μs <

 Memory estimate: 25.75 KiB, allocs estimate: 6.

10k
BenchmarkTools.Trial: 2464 samples with 1 evaluation.
 Range (min … max):  1.820 ms … 11.215 ms  ┊ GC (min … max): 0.00% … 81.53%
 Time  (median):     1.869 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.027 ms ±  1.191 ms  ┊ GC (mean ± σ):  7.49% ± 10.45%

  █▂                                                          
  ██▆▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ █
  1.82 ms      Histogram: log(frequency) by time     11.2 ms <

 Memory estimate: 2.91 MiB, allocs estimate: 38122.

 BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  46.833 μs …  3.458 ms  ┊ GC (min … max): 0.00% … 97.19%
 Time  (median):     48.458 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   50.596 μs ± 75.877 μs  ┊ GC (mean ± σ):  3.31% ±  2.17%

         ▂▃█▅▄▃ ▁                                              
  ▁▁▂▂▃▅█████████▆▅▄▄▄▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  46.8 μs         Histogram: frequency by time        54.7 μs <

 Memory estimate: 100.00 KiB, allocs estimate: 6.

 BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  27.667 μs …  3.601 ms  ┊ GC (min … max): 0.00% … 97.15%
 Time  (median):     29.042 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   33.599 μs ± 69.904 μs  ┊ GC (mean ± σ):  4.55% ±  2.19%

  ▃▇█▇▅▃▃▃▁               ▁▁                   ▁▁▁▁▁ ▁ ▁▁     ▂
  █████████▇▆▇▇▇▅▆▆▇▇▇▇▆████████▅▅▅▄▅▅▆▆▆▄▄▄▇████████████████ █
  27.7 μs      Histogram: log(frequency) by time        55 μs <

 Memory estimate: 100.00 KiB, allocs estimate: 6.

70k
BenchmarkTools.Trial: 353 samples with 1 evaluation.
 Range (min … max):  12.757 ms … 32.898 ms  ┊ GC (min … max): 0.00% … 58.64%
 Time  (median):     13.106 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   14.217 ms ±  4.491 ms  ┊ GC (mean ± σ):  7.59% ± 13.57%

  ██▂                                                       ▁  
  ███▄▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆█ ▆
  12.8 ms      Histogram: log(frequency) by time      32.7 ms <

 Memory estimate: 20.19 MiB, allocs estimate: 264625.

 BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  320.167 μs …   9.202 ms  ┊ GC (min … max): 0.00% … 95.58%
 Time  (median):     326.167 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   345.778 μs ± 172.046 μs  ┊ GC (mean ± σ):  1.90% ±  3.79%

  ▅█▇▅▄▂▂▁                                            ▁▁        ▁
  █████████▇▇▇▇▆▅▅▅▄▄▄▄▄▄▄▆▇█▇█▇▇▅▅▅▅▂▄▃▃▂▂▄▄▄▂▄▆▆█▇▇█████▇▆▆▇▆ █
  320 μs        Histogram: log(frequency) by time        496 μs <

 Memory estimate: 689.81 KiB, allocs estimate: 6.

 BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  190.500 μs …   3.536 ms  ┊ GC (min … max): 0.00% … 87.92%
 Time  (median):     196.917 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   206.381 μs ± 133.642 μs  ┊ GC (mean ± σ):  2.65% ±  3.91%

   ▇█▆▄▃▁ ▁                                                     ▂
  ██████████▇▇▅▃▄▄▃▄▄▃▃▃▁▄▃▁▁▁▁▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▆▆▇▆▇▆▆▇ █
  190 μs        Histogram: log(frequency) by time        347 μs <

 Memory estimate: 689.81 KiB, allocs estimate: 6.

WINNER: my_mul_mt_1! and my_mul_mt_2! (basically the same functions)

"""

##############################################################################
# CASE 2: PTDF and ABA MATRICES ##############################################

data = PowerFlowData(PTDFDCPowerFlow(), sys)
data_1 = deepcopy(data)
data_2 = deepcopy(data)

power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
matrix_data = data.power_network_matrix.data

@benchmark PF.my_mul_mt!(data.branch_flow_values, matrix_data, power_injection)
@benchmark PF.my_mul_mt_1!(data_1.branch_flow_values, matrix_data, power_injection)
@benchmark PF.my_mul_mt_2!(data_2.branch_flow_values, matrix_data, power_injection)

PF.my_mul_mt!(data.branch_flow_values, matrix_data, power_injection)
PF.my_mul_mt_1!(data_1.branch_flow_values, matrix_data, power_injection)
PF.my_mul_mt_2!(data_2.branch_flow_values, matrix_data, power_injection)

isapprox(data.branch_flow_values, data_2.branch_flow_values, atol = 1e-5)
isapprox(data.branch_flow_values, data_1.branch_flow_values, atol = 1e-6)
isapprox(data_1.branch_flow_values, data_2.branch_flow_values, atol = 1e-6)
isapprox(data_1.branch_flow_values, data.branch_flow_values, atol = 1e-6)

sum(abs.(data.branch_flow_values - data_2.branch_flow_values))
sum(abs.(data.branch_flow_values - data_1.branch_flow_values))
sum(abs.(data_2.branch_flow_values - data_1.branch_flow_values))

"""
2k
BenchmarkTools.Trial: 350 samples with 1 evaluation.
 Range (min … max):   6.203 ms … 26.308 ms  ┊ GC (min … max):  0.00% … 41.79%
 Time  (median):     15.606 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   14.296 ms ±  5.526 ms  ┊ GC (mean ± σ):  17.18% ± 21.12%

  ▆▁▂              ▁    ▁ ▁▂ ▃█▂     ▁                         
  █████▅▅██▅▄▅█▁▁▆██▇▆▆▁█▇██▄███▄▄▁▁▅█▆▅▅▇▄▅▁▅▇█▇▆█▅▄▅▅▇▇▇▅▅▆ ▆
  6.2 ms       Histogram: log(frequency) by time        26 ms <

 Memory estimate: 49.31 MiB, allocs estimate: 3207.

 BenchmarkTools.Trial: 864 samples with 1 evaluation.
 Range (min … max):  5.767 ms … 5.839 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     5.777 ms             ┊ GC (median):    0.00%
 Time  (mean ± σ):   5.778 ms ± 9.854 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

   ▅▂    █▅                                                  
  ▃██▃▂▂▅██▃▂▂▄▅▇▅▃▄▆█▆▄▃▄▃▃▃▄▄▂▃▃▃▃▃▄▃▃▃▃▃▂▂▃▃▂▂▃▂▂▂▂▁▂▂▁▂ ▃
  5.77 ms        Histogram: frequency by time       5.81 ms <

 Memory estimate: 25.34 KiB, allocs estimate: 3.

 BenchmarkTools.Trial: 8292 samples with 1 evaluation.
 Range (min … max):  504.208 μs …   6.224 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     528.208 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   592.610 μs ± 199.911 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▆█▄▃▃▃▃▃▂▂▂                                                   ▁
  ████████████▇▇▆▆▆▆▅▅▆▆▆▆▅▅▅▆▅▄▅▅▅▆▅▆▆▆▆▆▆▇▇▆▆▆▄▅▄▅▄▅▄▄▄▅▄▄▄▃▂ █
  504 μs        Histogram: log(frequency) by time       1.52 ms <

 Memory estimate: 25.34 KiB, allocs estimate: 3.

10k
BenchmarkTools.Trial: 21 samples with 1 evaluation.
 Range (min … max):  145.487 ms … 322.271 ms  ┊ GC (min … max):  0.00% … 11.74%
 Time  (median):     239.008 ms               ┊ GC (median):    14.15%
 Time  (mean ± σ):   238.210 ms ±  64.406 ms  ┊ GC (mean ± σ):   8.99% ±  6.84%

  ▃                                                    █    ▃    
  █▇▇▁▁▇▁▁▁▁▁▁▁▁▇▁▁▁▁▇▇▁▁▁▁▇▇▁▁▁▁▁▇▁▇▁▁▁▁▁▁▁▁▁▁▁▁▇▇▁▁▁▁█▁▁▇▁█▁▇ ▁
  145 ms           Histogram: frequency by time          322 ms <

 Memory estimate: 969.97 MiB, allocs estimate: 25413.

 BenchmarkTools.Trial: 43 samples with 1 evaluation.
 Range (min … max):  117.694 ms … 120.529 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     118.031 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   118.224 ms ± 497.927 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    █▆▁▃      ▁ ▁      ▃                                         
  ▄▁████▇▄▁▄▁▁█▄█▇▄▇▄▁▁█▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄ ▁
  118 ms           Histogram: frequency by time          121 ms <

 Memory estimate: 99.59 KiB, allocs estimate: 3.

 BenchmarkTools.Trial: 473 samples with 1 evaluation.
 Range (min … max):   9.640 ms …  13.237 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     10.545 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   10.546 ms ± 362.557 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

                          ▁▃▅▆▄█▇█▃▃ ▁                          
  ▂▃▂▃▂▂▂▃▄▄▅▃▃▄▃▆▅▅▅▆▆▇▆███████████▇██▄▆▄▆▅▃▄▂▄▁▂▃▂▂▃▂▁▂▂▂▄▂▂ ▄
  9.64 ms         Histogram: frequency by time         11.5 ms <

 Memory estimate: 99.59 KiB, allocs estimate: 3.

70k

"""

valid_ix = setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
p_inj = power_injection[valid_ix]
data.bus_angle[valid_ix] = data.aux_network_matrix.K \ p_inj