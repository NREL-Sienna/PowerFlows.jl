struct TxSteppingResidual
    data::PowerFlowData
    state::TxSteppingSolverState
    Rv::Vector{Float64}
    ref_mag::Vector{Float64}
end

function compute_residual!(
    residual::TxSteppingResidual,
    time_step::Int,
    lambda::Float64,
)
    data = residual.data
    state = residual.state
    V = state.V
    y_lambda = data.power_network_matrix.y_lambda
    network_I = y_lambda * V

    n = size(data.bus_type, 1)
    bus_types = @view data.bus_type[:, time_step]
    pv_mask = bus_types .== (PSY.ACBusTypes.PV,)

    # Power injections scaled by (1 - lambda): zero at lambda=1, full at lambda=0.
    S =
        (1 - lambda) .* (
            data.bus_active_power_injections[:, time_step] .-
            data.bus_active_power_withdrawals[:, time_step] .+
            im .* (
                data.bus_reactive_power_injections[:, time_step] .-
                data.bus_reactive_power_withdrawals[:, time_step]
            )
        )
    # handle Q_g's at PV buses. assumes that there's 0's in q_g for non-PVs.
    @assert all(state.q_g[.!pv_mask] .== 0.0) "Nonzero q_g values for non-PV buses not supported"
    S_imag = @view reinterpret(Float64, S)[2:2:end]
    S_imag .+= state.q_g

    device_I = conj.(S ./ V)
    I_net = network_I .- device_I
    residuals_Is = reinterpret(ComplexF64, @view(residual.Rv[1:(2 * n)]))
    residuals_Is .= I_net

    # Voltage constraints: one equation per PV bus.
    # Vset interpolated: (1 - lambda) * Vset_original + lambda * V_ref_mag
    Vset_mag =
        (1 - lambda) .* data.bus_magnitude[pv_mask, time_step] .+
        lambda .* residual.ref_mag[pv_mask]
    residual.Rv[(2 * n + 1):end] .= Vset_mag .^ 2 .- abs2.(V[pv_mask])

    # REF bus constraints: voltage fixed at setpoint
    for (ix, bt) in enumerate(bus_types)
        bt != PSY.ACBusTypes.REF && continue
        V_set = data.bus_magnitude[ix, time_step] * exp(im * data.bus_angles[ix, time_step])
        residual.Rv[2 * ix - 1] = real(state.V[ix]) - real(V_set)
        residual.Rv[2 * ix] = imag(state.V[ix]) - imag(V_set)
    end
end
