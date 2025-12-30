"""
    penalty_factors_brute_force(data::PowerFlowData; kwargs...)

Calculate the penalty factors for each bus in the power flow data using a brute force method.
This function calculates the penalty factors for each bus in the power flow data by perturbing the active power injection 
at each bus by a small step size and solving the power flow equations. 
The loss factor value is computed as the change in the reference bus power injection divided by the step size.

# Arguments
- `data::PowerFlowData`: The power flow data containing bus types, active power injections, and other relevant information.
- `step_size::Float64 = 1e-6`: The step size used to perturb the active power injection at each bus.
- `kwargs...`: Additional keyword arguments to be passed to the `solve_power_flow!` function.

# Returns
- `loss_factors::Array{Float64, 2}`: A 2D array of penalty factors for each bus and time step.

# Notes
- The reference bus type is assumed to remain constant between time steps.
- The initial power flow solution is computed to establish the starting power injection at the slack bus.
"""

function penalty_factors_brute_force(
    data::PowerFlowData,
    pf::ACPowerFlow;
    step_size::Float64 = 1e-6,
    kwargs...,
)
    if data.calculate_loss_factors
        @warn "Data with `calculate_loss_factors = true` passed to `penalty_factors_brute_force`:" *
              " this will re-compute the loss factors repeatedly, for no reason."
    end
    # we assume that the bus type for ref bus does not change between time steps
    ref, = PowerFlows.bus_type_idx(data, 1, (PSY.ACBusTypes.REF,))

    n_buses = first(size(data.bus_type))
    time_steps = collect(values(data.time_step_map))

    loss_factors = zeros(Float64, n_buses, length(time_steps))

    # initial PF to establish the ref power value
    solve_power_flow!(data; pf = pf, kwargs...)

    ref_power = copy(sum(data.bus_active_power_injections[ref, :]; dims = 1))

    for bx in 1:n_buses
        if bx in ref
            loss_factors[bx, :] .= 1.0
            continue
        end
        data.bus_active_power_injections[bx, :] .+= step_size
        solve_power_flow!(data; pf = pf, kwargs...)
        loss_factors[[bx], :] .=
            (sum(data.bus_active_power_injections[ref, :]; dims = 1) .- ref_power) ./
            step_size
        data.bus_active_power_injections[bx, :] .-= step_size
    end
    return loss_factors
end
