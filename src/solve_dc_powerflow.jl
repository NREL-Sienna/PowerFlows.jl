"""
    solve_powerflow!(data::PTDFPowerFlowData)

Evaluates the PTDF power flow and writes the result to the fields of the 
[`PTDFPowerFlowData`](@ref) structure.

This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_activepower_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_activepower_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
function solve_powerflow!(
    data::PTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    # get net power injections
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    # evaluate flows
    data.branch_activepower_flow_from_to .=
        data.power_network_matrix.data' * power_injection
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
    # evaluate bus angles
    p_inj = power_injection[data.valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, :] .= p_inj
    data.converged .= true
    return
end

"""
    solve_powerflow!(data::vPTDFPowerFlowData)

Evaluates the virtual PTDF power flow and writes the results to the fields 
of the [`vPTDFPowerFlowData`](@ref) structure.


This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_activepower_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_activepower_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
function solve_powerflow!(
    data::vPTDFPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.aux_network_matrix.data)
    full_factor!(solver_cache, data.aux_network_matrix.data)
    power_injection = data.bus_activepower_injection .- data.bus_activepower_withdrawals
    data.branch_activepower_flow_from_to .=
        my_mul_mt(data.power_network_matrix, power_injection)
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
    p_inj = power_injection[data.valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, :] .= p_inj
    data.converged .= true
    return
end

# TODO: solve just for some lines with vPTDF

"""
    solve_powerflow!(data::ABAPowerFlowData)

Evaluates the DC power flow and writes the results (branch flows) to the fields 
of the [`ABAPowerFlowData`](@ref) structure.


This function modifies the following fields of `data`, setting them to the computed values:
- `data.bus_angles`: the bus angles for each bus in the system.
- `data.branch_activepower_flow_from_to`: the active power flow from the "from" bus to the "to" bus of each branch
- `data.branch_activepower_flow_to_from`: the active power flow from the "to" bus to the "from" bus of each branch

Additionally, it sets `data.converged` to `true`, indicating that the power flow calculation was successful.
"""
# DC flow: ABA and BA case
function solve_powerflow!(
    data::ABAPowerFlowData,
)
    solver_cache = KLULinSolveCache(data.power_network_matrix.data)
    full_factor!(solver_cache, data.power_network_matrix.data)
    # get net injections
    power_injection = data.bus_activepower_injection - data.bus_activepower_withdrawals
    # save angles and power flows
    p_inj = power_injection[data.valid_ix, :]
    solve!(solver_cache, p_inj)
    data.bus_angles[data.valid_ix, :] .= p_inj
    data.branch_activepower_flow_from_to .= data.aux_network_matrix.data' * data.bus_angles
    data.branch_activepower_flow_to_from .= -data.branch_activepower_flow_from_to
    data.converged .= true
    return
end

# SINGLE PERIOD ##############################################################

"""
    solve_powerflow(
        ::T,
        sys::PSY.System;
    ) where T <: Union{PTDFDCPowerFlow, vPTDFDCPowerFlow, DCPowerFlow}


Evaluates the provided DC power flow method `T` on the `system`, returning a dictionary of 
`DataFrame`s containing the calculated branch flows and bus angles.

Provided for convenience: this interface bypasses the need to create a `PowerFlowData` 
struct, but that's still what's happening under the hood.

# Example
```julia
using PowerFlows, PowerSystemCaseBuilder
sys = PowerSystemCaseBuilder.build_system(PSB.PSITestSystems, "c_sys5")
d = solve_powerflow(DCPowerFlow(), sys)
display(d["1"]["flow_results"])
display(d["1"]["bus_results"])
```
"""
function solve_powerflow(
    ::T,
    sys::PSY.System;
) where {T <: Union{PTDFDCPowerFlow, vPTDFDCPowerFlow, DCPowerFlow}}
    data = PowerFlowData(T(), sys)
    solve_powerflow!(data)
    return write_results(data, sys)
end

# MULTI PERIOD ###############################################################

"""
    solve_powerflow(
        data::T,
        sys::PSY.System;
    ) where T <: Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData}

Evaluates the DC power flow associated with `T`, writes the result to `data`, and returns 
a dictionary of `DataFrame`s containing the calculated branch flows and bus angles. 

Note that `data` must have been created from the system `sys` using one of the 
[`PowerFlowData`](@ref) constructors, such as `PTDFPowerFlowData(sys)` or 
`vPTDFPowerFlowData(sys)`.

# Example
```julia
using PowerFlows, PowerSystemCaseBuilder
sys = PowerSystemCaseBuilder.build_system(PSITestSystems, "c_sys14")
data = PowerFlowData(PTDFDCPowerFlow(), sys, time_steps = 2)
d = solve_powerflow(data, sys)
display(d["2"]["flow_results"])
```
"""
function solve_powerflow(
    data::T,
    sys::PSY.System;
) where {T <: Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData}}
    solve_powerflow!(data)
    return write_results(data, sys)
end
