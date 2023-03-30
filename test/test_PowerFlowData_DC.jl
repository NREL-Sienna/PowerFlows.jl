## get the needed packges
using Revise
using Test
using Logging
using PowerSystems
using PowerSystemCaseBuilder
using InfrastructureSystems
using PowerFlows
using LinearAlgebra
using CSV
using DataFrames
using PowerNetworkMatrices
import KLU: klu
const PNM = PowerNetworkMatrices
const IS = InfrastructureSystems
const PSB = PowerSystemCaseBuilder
const PSY = PowerSystems

## load the sytem (RTS) ######################################################
sys = PSB.build_system(PSITestSystems, "test_RTS_GMLC_sys")

# get the AC branches
branches = PNM.get_ac_branches(sys)
buses = PNM.get_buses(sys)

# get the branch names and bus numbers
line_ax = [PSY.get_name(branch) for branch in branches]
bus_ax = [PSY.get_number(bus) for bus in buses]

# create dictionaries for reference
bus_ix = PNM.make_ax_ref(bus_ax)
branch_ix = PNM.make_ax_ref(line_ax)

# get injections' active and reactive power
# ? just sum injections on the same buses
# ? how to treat buses with multiple injections
bus_number_injections = [PSY.get_number(PSY.get_bus(d)) for d in PSY.get_components(StaticInjection, sys)]
bus_idx = sortperm(bus_number_injections)
bus_activepower_injection = [PSY.get_active_power(d) for d in PSY.get_components(StaticInjection, sys)][bus_idx]
bus_reactivepower_injection = [PSY.get_reactive_power(d) for d in PSY.get_components(StaticInjection, sys)][bus_idx]

# get bus type
bus_type = [PSY.get_bustype(b) for b in buses]

# bus angle
bus_angle = [PSY.get_angle(b) for b in PSY.get_components(PSY.Bus, sys)]

# bus magnitude
bus_magnitude = [PSY.get_magnitude(b) for b in PSY.get_components(PSY.Bus, sys)]

## evaluate the different matrices ###########################################

# A
A = PNM.IncidenceMatrix(sys)

# BA
BA = PNM.BA_Matrix(sys)

# ABA # ! not working at the moment, wait for PNM to be updated
ABA = A.data[:, setdiff(1:end, A.ref_bus_positions)]' * BA.data
# ABA = PNM.ABA_Matrix(sys)

# LU factorization matrices for ABA
K = klu(ABA)

# get vector of angles theta (inv(ABA)*bus_activepower_injection)
theta = K \ bus_activepower_injection

# PTDF
ptdf = PNM.PTDF(sys)

## compute power flow with different methods #################################

# compute power flow using PTDF matrix

# compute power flow using B and theta

## store data in PowerFlowData struc
