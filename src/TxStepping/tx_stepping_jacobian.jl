struct TxSteppingJacobian
    data::PowerFlowData
    Jv::SparseArrays.SparseMatrixCSC{Float64, J_INDEX_TYPE}
    rhs::Vector{Float64}
    state::TxSteppingSolverState
    ref_mag::Vector{Float64}
end

"""
Build the sparse Jacobian matrix with the correct sparsity pattern for the TxStepping
formulation. The Jacobian has dimension `(2n + n_pv) x (2n + n_pv)`:
- Upper-left `2n x 2n` block mirrors y_lambda's structure (each complex entry → 2x2 real block)
- Extra rows/columns for PV bus voltage magnitude constraints
"""
function _build_tx_stepping_jacobian(
    data::PowerFlowData,
    state::TxSteppingSolverState,
    ref_mag::Vector{Float64},
    time_step::Int = 1,
)
    y_lambda = data.power_network_matrix.y_lambda
    n = size(data.bus_type, 1)
    bus_types = @view data.bus_type[:, time_step]
    pv_indices =
        J_INDEX_TYPE[ix for (ix, bt) in enumerate(bus_types) if bt == PSY.ACBusTypes.PV]
    n_pv = length(pv_indices)
    dim = 2 * n + n_pv

    # Build the 2n x 2n block from y_lambda's sparsity pattern.
    # Each (i, j) entry in y_lambda → four entries at (2i-1, 2j-1), (2i-1, 2j), (2i, 2j-1), (2i, 2j).
    y_I, y_J, _ = SparseArrays.findnz(y_lambda)
    nnz_y = length(y_I)

    # Total entries: 4 per y_lambda nonzero + 4 per PV bus
    n_total = 4 * nnz_y + 4 * n_pv
    all_I = Vector{J_INDEX_TYPE}(undef, n_total)
    all_J = Vector{J_INDEX_TYPE}(undef, n_total)

    # Y-block: vectorized index expansion.
    # For each nonzero (y_I[k], y_J[k]), emit a 2x2 block at rows (2i-1, 2i) × cols (2j-1, 2j).
    y_rows = J_INDEX_TYPE.(y_I)
    y_cols = J_INDEX_TYPE.(y_J)
    yblock = @view all_I[1:(4 * nnz_y)]
    jblock = @view all_J[1:(4 * nnz_y)]
    # Lay out as: [odd rows; even rows; odd rows; even rows] × [odd cols; odd cols; even cols; even cols]
    yblock[1:4:end] .= 2 .* y_rows .- 1  # (2i-1, 2j-1)
    yblock[2:4:end] .= 2 .* y_rows        # (2i,   2j-1)
    yblock[3:4:end] .= 2 .* y_rows .- 1  # (2i-1, 2j)
    yblock[4:4:end] .= 2 .* y_rows        # (2i,   2j)
    jblock[1:4:end] .= 2 .* y_cols .- 1
    jblock[2:4:end] .= 2 .* y_cols .- 1
    jblock[3:4:end] .= 2 .* y_cols
    jblock[4:4:end] .= 2 .* y_cols

    # PV bus coupling entries: for each PV bus ix with q_g column q_col,
    # add 4 entries: (2ix-1, q_col), (2ix, q_col), (q_col, 2ix-1), (q_col, 2ix).
    pv_offset = 4 * nnz_y
    q_cols = J_INDEX_TYPE.(2 * n .+ (1:n_pv))
    bus_2i_minus1 = J_INDEX_TYPE.(2 .* pv_indices .- 1)
    bus_2i = J_INDEX_TYPE.(2 .* pv_indices)
    pv_I = @view all_I[(pv_offset + 1):end]
    pv_J = @view all_J[(pv_offset + 1):end]
    pv_I[1:4:end] .= bus_2i_minus1
    pv_J[1:4:end] .= q_cols   # (2i-1, q_col)
    pv_I[2:4:end] .= bus_2i
    pv_J[2:4:end] .= q_cols   # (2i,   q_col)
    pv_I[3:4:end] .= q_cols
    pv_J[3:4:end] .= bus_2i_minus1  # (q_col, 2i-1)
    pv_I[4:4:end] .= q_cols
    pv_J[4:4:end] .= bus_2i         # (q_col, 2i)

    all_V = zeros(Float64, n_total)
    Jv = SparseArrays.sparse(all_I, all_J, all_V, dim, dim)
    rhs = zeros(Float64, dim)
    return TxSteppingJacobian(data, Jv, rhs, state, ref_mag)
end

function stamp_jacobian!(
    Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE},
    rhs::Vector{<:Number},
    V::Complex,
    S::Complex,
    i::Int,
    is_pv::Bool,
    q_index::Int,
    Vset::Float64,
    q_g::Float64,
)
    # Device current I_dev = conj(S/V) and its Jacobian w.r.t. [vr, vi].
    # g = S/V² encodes the 2x2 Jacobian block: diagonal ±Re(g), off-diagonal Im(g).
    I_dev = conj(S / V)
    g = S / V^2

    Jv[2 * i - 1, 2 * i - 1] += real(g)
    Jv[2 * i - 1, 2 * i] -= imag(g)
    Jv[2 * i, 2 * i - 1] -= imag(g)
    Jv[2 * i, 2 * i] -= real(g)

    # RHS = I_dev + J_device * V evaluates to 2*I_dev because J_device * V = I_dev.
    rhs_val = 2 * I_dev

    if is_pv
        dIdQ = conj(im / V)

        Jv[2 * i - 1, q_index] = -real(dIdQ)
        Jv[2 * i, q_index] = -imag(dIdQ)
        Jv[q_index, 2 * i - 1] = -2 * real(V)
        Jv[q_index, 2 * i] = -2 * imag(V)

        rhs_val -= dIdQ * q_g
        rhs[q_index] -= Vset^2 + abs2(V)
    end

    rhs[2 * i - 1] += real(rhs_val)
    rhs[2 * i] += imag(rhs_val)
end

function _update_jacobian!(
    Jv::SparseArrays.SparseMatrixCSC{Float64, J_INDEX_TYPE},
    rhs::Vector{<:Number},
    state::TxSteppingSolverState,
    data::PowerFlowData,
    ref_mag::Vector{Float64},
    time_step::Int,
    lambda::Float64,
)
    # write network backbone to Jacobian. assumptions:
    # (1) y lambda has already been updated
    # (2) upper-left 2n-by-2n block of Jv has same sparse structure as y_lambda,
    #     except with a 2x2 block of entries for each complex entry in y_lambda.
    y_lambda = data.power_network_matrix.y_lambda
    n = size(data.bus_type, 1)
    y_lambda_nonzeros = SparseArrays.nonzeros(y_lambda)
    Jv_nonzeros = SparseArrays.nonzeros(Jv)

    for i in 1:n
        y_col_range = SparseArrays.nzrange(y_lambda, i)
        y_col_len = length(y_col_range)

        # Reinterpret ComplexF64 column as Float64 pairs [g₁, b₁, g₂, b₂, ...].
        # Odd Jv column (2i-1) gets [g, b] per entry; even column (2i) gets [-b, g].
        y_real = reinterpret(Float64, @view(y_lambda_nonzeros[y_col_range]))
        # only write first (2*y_col_len-1) items in column: skip the PV coupling entries
        odd_dest =
            @view Jv_nonzeros[Jv.colptr[2 * i - 1]:(Jv.colptr[2 * i - 1] + 2 * y_col_len - 1)]
        even_dest =
            @view Jv_nonzeros[Jv.colptr[2 * i]:(Jv.colptr[2 * i] + 2 * y_col_len - 1)]
        copyto!(odd_dest, y_real)                            # [g₁, b₁, g₂, b₂, ...]
        copyto!(even_dest, y_real)
        reinterpret(ComplexF64, even_dest) .*= im            # [g, b] → [-b, g]
    end
    rhs .= 0.0

    # write power injection contributions to Jacobian and update rhs.
    q_index = 2 * n + 1
    bus_types = @view data.bus_type[:, time_step]
    for (ix, bt) in enumerate(bus_types)
        bt == PSY.ACBusTypes.REF && continue
        S_original =
            (1 - lambda) * (
                data.bus_active_power_injections[ix, time_step] -
                data.bus_active_power_withdrawals[ix, time_step] +
                im * (
                    data.bus_reactive_power_injections[ix, time_step] -
                    data.bus_reactive_power_withdrawals[ix, time_step]
                )
            )
        S = S_original + state.q_g[ix] * im
        is_pv = bt == PSY.ACBusTypes.PV
        Vset = if is_pv
            (1 - lambda) * data.bus_magnitude[ix, time_step] + lambda * ref_mag[ix]
        else
            0.0
        end
        stamp_jacobian!(Jv, rhs, state.V[ix], S, ix, is_pv, q_index, Vset, state.q_g[ix])
        is_pv && (q_index += 1)
    end

    # REF bus rows: overwrite with identity (V_ref = V_set).
    # Zero all entries in REF bus rows (the Y copy wrote admittances there),
    # then set the diagonal to 1.
    # PERF: store ahead of time instead of re-calculating each time.
    ref_rows = Set{J_INDEX_TYPE}()
    for (ix, bt) in enumerate(bus_types)
        bt != PSY.ACBusTypes.REF && continue
        push!(ref_rows, J_INDEX_TYPE(2 * ix - 1))
        push!(ref_rows, J_INDEX_TYPE(2 * ix))
    end
    Jv_rowvals = SparseArrays.rowvals(Jv)
    for k in eachindex(Jv_rowvals)
        if Jv_rowvals[k] in ref_rows
            Jv_nonzeros[k] = 0.0
        end
    end
    for (ix, bt) in enumerate(bus_types)
        bt != PSY.ACBusTypes.REF && continue
        Jv[2 * ix - 1, 2 * ix - 1] = 1.0
        Jv[2 * ix, 2 * ix] = 1.0
        V_set = data.bus_magnitude[ix, time_step] *
                exp(im * data.bus_angles[ix, time_step])
        rhs[2 * ix - 1] += real(V_set)
        rhs[2 * ix] += imag(V_set)
    end
end
