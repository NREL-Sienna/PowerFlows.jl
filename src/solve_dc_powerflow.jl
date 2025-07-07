const PTDFPowerFlowData = PowerFlowData{
    PNM.PTDF{
        Tuple{Vector{Int64}, Vector{String}},
        Tuple{Dict{Int64, Int64}, Dict{String, Int64}},
        Matrix{Float64},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

const vPTDFPowerFlowData = PowerFlowData{
    PNM.VirtualPTDF{
        Tuple{Vector{String}, Vector{Int64}},
        Tuple{Dict{String, Int64}, Dict{Int64, Int64}},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

const ABAPowerFlowData = PowerFlowData{
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
    PNM.BA_Matrix{
        Tuple{Vector{Int64}, Vector{String}},
        Tuple{Dict{Int64, Int64}, Dict{String, Int64}}},
}

"""
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::PTDFPowerFlowData`:
        PTDFPowerFlowData structure containing all the information related to the system's power flow.
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
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::vPTDFPowerFlowData`:
        vPTDFPowerFlowData structure containing all the information related to the system's power flow.
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
Evaluates the power flows on each system's branch and updates the PowerFlowData structure.

# Arguments:
- `data::ABAPowerFlowData`:
        ABAPowerFlowData structure containing all the information related to the system's power flow.
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
Evaluates the power flows on the system's branches by means of the PTDF, virtual PTDF,
or DC power flow method: the type first parameter (a `PTDFDCPowerFlow`, `vPTDFDCPowerFlow`, 
or `DCPowerFlow`) selects the method to be used. Returns a dictionary containing a 
`DataFrame` for the single timestep considered, storing the branch flows and bus 
voltages for the input `PSY.System`.

# Arguments:
- `::Union{PTDFDCPowerFlow, vPTDFDCPowerFlow, DCPowerFlow}`:
        the method of power flow evaluation to be used.
- `sys::PSY.System`:
        container gathering the system data used for the evaluation of flows
        and angles.
"""
function solve_powerflow(
    ::T,
    sys::PSY.System;
    correct_bustypes::Bool = false,
) where {T <: AbstractDCPowerFlow}
    data = PowerFlowData(T(), sys; correct_bustypes = correct_bustypes)
    solve_powerflow!(data)
    return write_results(data, sys)
end

# MULTI PERIOD ###############################################################

"""
Evaluates the power flows on the system's branches by means of the method associated with
the `PowerFlowData` structure `data`, which can be one of `PTDFPowerFlowData`,
`vPTDFPowerFlowData`, or `ABAPowerFlowData`.
Returns a dictionary of `DataFrame`s, each containing the branch flows and bus voltages for
the input `PSY.System` at that timestep.

# Arguments:
- `data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData}`:
        `PowerFlowData` structure containing the system's data per each timestep
        considered, as well as the associated matrix for the power flow.
- `sys::PSY.System`:
        container gathering the system data.
"""
function solve_powerflow(
    data::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData},
    sys::PSY.System;
)
    solve_powerflow!(data)
    return write_results(data, sys)
end
