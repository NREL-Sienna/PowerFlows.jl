"""
DC contingency analysis — the Sienna counterpart of PSS/e DCCC.

For each contingency (an AC-branch outage described by a `PSY.Outage`
supplemental attribute or a pre-built `PNM.NetworkModification`), this module
computes the post-contingency bus angles and arc flows across all time steps,
reusing linear algebra work from the base factorization.

Strategies (Phase 1 = refactor only):
- `:refactor` — numeric refactor of `ABA_base + ΔABA`, shared symbolic
  factorization from the base.
- `:woodbury` — Sherman-Morrison-Woodbury correction via `ABA_base⁻¹ · U · W⁻¹ · Uᵀ`
  (Phase 2).
- `:auto` — Woodbury if `M = length(mod.arc_modifications) ≤ woodbury_max_rank`,
  else refactor (Phase 2).

Timesteps are batched as a single multi-RHS KLU solve per contingency.
"""

# ---------------------------------------------------------------------------
# Contingency → NetworkModification resolution
# ---------------------------------------------------------------------------

"""
    _resolve_modifications(aba, sys, contingencies) -> (mods, labels)

Accept a vector of `PSY.Outage` or `PNM.NetworkModification` and return the
aligned vectors of modifications and human-readable labels. `ContingencySpec`
wrapping is internal and not exposed.
"""
function _resolve_modifications(
    ba::PNM.BA_Matrix,
    sys::PSY.System,
    contingencies::AbstractVector{<:PSY.Outage},
)
    mods = Vector{PNM.NetworkModification}(undef, length(contingencies))
    labels = Vector{String}(undef, length(contingencies))
    for (i, outage) in enumerate(contingencies)
        # NetworkModification uses the network-reduction data on `ba` (shared
        # with ABA) plus `ba`'s arc_lookup to resolve outage components → arcs.
        mods[i] = PNM.NetworkModification(ba, sys, outage)
        labels[i] = _contingency_label(outage, i)
    end
    return mods, labels
end

function _resolve_modifications(
    ::PNM.BA_Matrix,
    ::PSY.System,
    contingencies::AbstractVector{<:PNM.NetworkModification},
)
    mods = collect(contingencies)
    labels = [isempty(m.label) ? "ctg_$(i)" : m.label for (i, m) in enumerate(mods)]
    return mods, labels
end

function _contingency_label(outage::PSY.Outage, i::Int)
    # PSY.Outage does not carry a user-facing `name`; use UUID for a stable,
    # unique label. Callers that want friendly labels should pass
    # `Vector{PNM.NetworkModification}` with the `label` field pre-populated.
    return "ctg_$(i)_$(string(IS.get_uuid(outage)))"
end

# ---------------------------------------------------------------------------
# Arc → (fb_ix, tb_ix) precomputation on the FULL bus axis
# ---------------------------------------------------------------------------

"""
    _arc_endpoint_indices(base_data) -> (fb_full, tb_full)

Return vectors of length `n_arcs` mapping each arc to its from/to bus positions
on the FULL bus axis (the axis of the result 3D arrays), not the reduced
(ref-bus-removed) axis. Used for flow reconstruction.
"""
function _arc_endpoint_indices(base_data::ABAPowerFlowData)
    arcs = get_arc_axis(base_data)
    bus_lookup = get_bus_lookup(base_data)
    n_arcs = length(arcs)
    fb = Vector{Int}(undef, n_arcs)
    tb = Vector{Int}(undef, n_arcs)
    @inbounds for (i, arc) in enumerate(arcs)
        fb[i] = bus_lookup[first(arc)]
        tb[i] = bus_lookup[last(arc)]
    end
    return fb, tb
end

# ---------------------------------------------------------------------------
# ABA nzval patching — in-place mutation to preserve sparsity pattern
# ---------------------------------------------------------------------------

"""
    _build_nzval_patch_map(ABA_base, aba_delta) -> Vector{Int}

For each `(row, col)` position mentioned in `aba_delta.delta`, return the index
into `ABA_base.nzval` of that entry. Used to mutate the base `nzval` in place
without allocating a new sparse matrix (preserves `colptr`/`rowval` so
`numeric_refactor!` with `check_pattern=true` stays valid).

Returns a vector of the same length as `SparseArrays.nnz(aba_delta.delta)` aligned with
iteration `(rows, cols, vals) = SparseArrays.findnz(aba_delta.delta)`.
"""
function _build_nzval_patch_map(
    ABA_base::SparseMatrixCSC{Float64, Int},
    aba_delta::PNM.ABADelta,
)
    delta_sp = aba_delta.delta
    rows = SparseArrays.rowvals(delta_sp)
    n_nz = SparseArrays.nnz(delta_sp)
    base_colptr = ABA_base.colptr
    base_rowval = ABA_base.rowval
    patch_idx = Vector{Int}(undef, n_nz)
    counter = 0
    for col in 1:size(delta_sp, 2)
        col_range = SparseArrays.nzrange(delta_sp, col)
        for nz in col_range
            counter += 1
            row = rows[nz]
            # Find this (row, col) in ABA_base's CSC structure
            base_range = base_colptr[col]:(base_colptr[col + 1] - 1)
            found = 0
            for bnz in base_range
                if base_rowval[bnz] == row
                    found = bnz
                    break
                end
            end
            if found == 0
                throw(
                    ArgumentError(
                        "ABA_base is missing the (row=$row, col=$col) structural entry " *
                        "needed by the delta. Modification introduces a new sparsity pattern " *
                        "position; cannot use in-place nzval patching.",
                    ),
                )
            end
            patch_idx[counter] = found
        end
    end
    return patch_idx
end

# ---------------------------------------------------------------------------
# Flow reconstruction with outaged arcs
# ---------------------------------------------------------------------------

"""
    _reconstruct_flows!(flow_out, BA_base, θ_full, mod_arc_indices, mod_delta_b,
        fb_full, tb_full, arc_susceptances)

Compute post-contingency arc flows into `flow_out` (n_arcs × n_ts).

For arcs NOT touched by the contingency: `flow = BA_baseᵀ · θ`.

For arcs in `mod_arc_indices`: the base BA row bakes in the pre-outage
susceptance; replace those rows with `b_eff · (θ_fb - θ_tb)` using the
post-contingency effective susceptance `b_eff = b_base + Δb`. For a full
outage (`Δb = -b_base`), `b_eff = 0` and flow is zero.
"""
function _reconstruct_flows!(
    flow_out::AbstractMatrix{Float64},
    BA_base::PNM.BA_Matrix,
    θ_full::AbstractMatrix{Float64},
    mod_arc_indices::Vector{Int},
    mod_delta_b::Vector{Float64},
    fb_full::Vector{Int},
    tb_full::Vector{Int},
    arc_susceptances::Vector{Float64},
)
    # flow_out := BA_baseᵀ · θ
    LinearAlgebra.mul!(flow_out, transpose(BA_base.data), θ_full)
    # Patch modified arcs with post-contingency effective susceptance
    @inbounds for (e_pos, arc_idx) in enumerate(mod_arc_indices)
        b_eff = arc_susceptances[arc_idx] + mod_delta_b[e_pos]
        if abs(b_eff) < eps()
            # full outage → zero flow
            for t in 1:size(flow_out, 2)
                flow_out[arc_idx, t] = 0.0
            end
        else
            ifrom = fb_full[arc_idx]
            ito = tb_full[arc_idx]
            for t in 1:size(flow_out, 2)
                flow_out[arc_idx, t] = b_eff * (θ_full[ifrom, t] - θ_full[ito, t])
            end
        end
    end
    return flow_out
end

# ---------------------------------------------------------------------------
# Base case → slab 1 population
# ---------------------------------------------------------------------------

"""
    _populate_base_slab!(result, base_data)

Copy base-case results into slab 1 AND broadcast the contingency-invariant
fields (bus magnitudes, bus types, injections, withdrawals) across ALL
contingency slabs at once. DC contingency analysis only changes branch flows
and angles — pumping base injections/types into every per-contingency
`_store_contingency_slab!` call would write hundreds of MB of identical data
on WECC-scale studies.
"""
function _populate_base_slab!(
    result::TimeContingencyPowerFlowData,
    base_data::ABAPowerFlowData,
)
    @inbounds begin
        # Slab 1 — base case (per-contingency varying fields)
        result.bus_angles[:, :, 1] .= base_data.bus_angles
        result.arc_active_power_flow_from_to[:, :, 1] .=
            base_data.arc_active_power_flow_from_to
        result.arc_active_power_flow_to_from[:, :, 1] .=
            base_data.arc_active_power_flow_to_from
        result.arc_angle_differences[:, :, 1] .= base_data.arc_angle_differences
        result.converged[:, 1] .= base_data.converged
        if result.arc_active_power_losses !== nothing &&
           base_data.arc_active_power_losses !== nothing
            result.arc_active_power_losses[:, :, 1] .= base_data.arc_active_power_losses
        end
        # Invariants — broadcast across ALL slabs (base + contingencies).
        # 3D `[:, :, k]` slab from a 2D source is `result_field[:, :, :] .= reshape(src, …, 1)`.
        n_ctg = size(result.bus_angles, 3)
        for k in 1:n_ctg
            result.bus_magnitude[:, :, k] .= base_data.bus_magnitude
            result.bus_type[:, :, k] .= base_data.bus_type
            result.bus_active_power_injections[:, :, k] .=
                base_data.bus_active_power_injections
            result.bus_reactive_power_injections[:, :, k] .=
                base_data.bus_reactive_power_injections
            result.bus_active_power_withdrawals[:, :, k] .=
                base_data.bus_active_power_withdrawals
            result.bus_reactive_power_withdrawals[:, :, k] .=
                base_data.bus_reactive_power_withdrawals
        end
    end
    return
end

# ---------------------------------------------------------------------------
# Angle diff helper (contingency slab)
# ---------------------------------------------------------------------------

function _compute_arc_angle_differences_slab!(
    angle_diff::AbstractMatrix{Float64},
    θ_full::AbstractMatrix{Float64},
    fb_full::Vector{Int},
    tb_full::Vector{Int},
)
    @inbounds for (i, (ifrom, ito)) in enumerate(zip(fb_full, tb_full))
        for t in 1:size(θ_full, 2)
            angle_diff[i, t] = θ_full[ifrom, t] - θ_full[ito, t]
        end
    end
    return
end

# ---------------------------------------------------------------------------
# Post-contingency connectivity: count electrical islands after removing
# the arcs touched by a modification. Cheap BFS over the base arc axis with
# the modified arcs filtered out (full outage) or scaled (partial).
# ---------------------------------------------------------------------------

"""
    _count_post_contingency_islands(base_data, mod) -> Int

Count the number of connected components in the post-contingency network.
Arcs in `mod.arc_modifications` are considered "removed" iff the outage is
full (`|Δb + b_base| < eps()`); partial outages leave the arc connected.
Returns `1` for a connected network, `>1` if the contingency creates islands.

O(n_bus + n_arcs) per call; called once per islanding contingency.
"""
function _count_post_contingency_islands(
    base_data::ABAPowerFlowData,
    mod::PNM.NetworkModification,
    arc_susceptances::Vector{Float64},
)
    ba = base_data.aux_network_matrix
    arc_ax = PNM.get_arc_axis(ba)
    bus_lookup = PNM.get_bus_lookup(ba)
    nr = PNM.get_network_reduction_data(ba)

    # Identify arcs "disconnected" by this contingency (full outage only).
    removed = falses(length(arc_ax))
    for am in mod.arc_modifications
        b_eff = arc_susceptances[am.arc_index] + am.delta_b
        if abs(b_eff) < eps()
            removed[am.arc_index] = true
        end
    end

    # Build an undirected adjacency list of remaining arcs.
    n_bus = length(PNM.get_bus_axis(ba))
    adj = [Int[] for _ in 1:n_bus]
    for (i, arc) in enumerate(arc_ax)
        removed[i] && continue
        fb = PNM.get_bus_index(arc[1], bus_lookup, nr)
        tb = PNM.get_bus_index(arc[2], bus_lookup, nr)
        push!(adj[fb], tb)
        push!(adj[tb], fb)
    end

    # BFS from each unvisited bus; count components.
    visited = falses(n_bus)
    n_components = 0
    queue = Int[]
    for start in 1:n_bus
        visited[start] && continue
        n_components += 1
        empty!(queue)
        push!(queue, start)
        visited[start] = true
        while !isempty(queue)
            u = popfirst!(queue)
            for v in adj[u]
                if !visited[v]
                    visited[v] = true
                    push!(queue, v)
                end
            end
        end
    end
    return n_components
end

# ---------------------------------------------------------------------------
# Per-task mutable state bundle (for parallelism via Channel pool)
# ---------------------------------------------------------------------------

"""
    ContingencyScratch

Bundle of per-worker mutable state for the contingency solve loop. One instance
per parallel task; distributed via a blocking `Channel` pool so that
`Threads.@threads :dynamic` can check out/in state without `threadid()`.

All fields are task-owned after `take!` — no field is shared across tasks.
"""
mutable struct ContingencyScratch
    # Per-task KLU cache for the refactor path (mutated every iteration).
    cache::KLULinSolveCache{Int64}
    # Per-task KLU factorization of the UNMODIFIED base ABA (never refactored).
    # Used by the Woodbury path so ldiv!(base_K, Z_buf) runs without racing on
    # the shared aba.K.common struct that klu_solve mutates.
    base_K::PNM.KLU.KLUFactorization{Float64, Int}
    nzval_scratch::Vector{Float64}
    ABA_patch_mat::SparseMatrixCSC{Float64, Int}
    rhs_buf::Matrix{Float64}
    θ_full_buf::Matrix{Float64}
    flow_buf::Matrix{Float64}
end

function ContingencyScratch(aba::PNM.ABA_Matrix, n_bus_full::Int, n_arcs::Int, n_ts::Int)
    nzval_base = copy(aba.data.nzval)
    # Patch matrix + refactor cache (refactor path)
    ABA_patch_mat = SparseMatrixCSC(
        size(aba.data, 1), size(aba.data, 2),
        copy(aba.data.colptr), copy(aba.data.rowval),
        copy(nzval_base),
    )
    cache = KLULinSolveCache(ABA_patch_mat)
    full_factor!(cache, ABA_patch_mat)
    # Independent base KLU factor for Woodbury path. A second copy of the ABA
    # sparse matrix is needed because KLUFactorization keeps references to the
    # passed colptr/rowval/nzval arrays.
    base_mat = SparseMatrixCSC(
        size(aba.data, 1), size(aba.data, 2),
        copy(aba.data.colptr), copy(aba.data.rowval),
        copy(nzval_base),
    )
    base_K = PNM.KLU.klu(base_mat)
    return ContingencyScratch(
        cache,
        base_K,
        similar(nzval_base),
        ABA_patch_mat,
        Matrix{Float64}(undef, size(aba.data, 1), n_ts),
        Matrix{Float64}(undef, n_bus_full, n_ts),
        Matrix{Float64}(undef, n_arcs, n_ts),
    )
end

# ---------------------------------------------------------------------------
# Per-contingency serial solve
# ---------------------------------------------------------------------------

"""
Solve one contingency (refactor path, Phase 1 MVP).

Arguments are all pre-allocated / shared across the loop to minimize allocation.
`cache` is mutated (new numeric factor); `nzval_scratch` is mutated in place.
"""
function _solve_one_contingency_refactor!(
    result::TimeContingencyPowerFlowData,
    k::Int,
    mod::PNM.NetworkModification,
    base_data::ABAPowerFlowData,
    aba::PNM.ABA_Matrix,
    scratch::ContingencyScratch,
    nzval_base::Vector{Float64},
    fb_full::Vector{Int},
    tb_full::Vector{Int},
    arc_susceptances::Vector{Float64},
    injections_full::Matrix{Float64},
    valid_ix::AbstractVector{Int},
)
    d = PNM.compute_aba_delta(aba, mod)
    copyto!(scratch.nzval_scratch, nzval_base)
    patch_map = _build_nzval_patch_map(aba.data, d)
    delta_vals = SparseArrays.nonzeros(d.delta)
    @inbounds for (i, pos) in enumerate(patch_map)
        scratch.nzval_scratch[pos] += delta_vals[i]
    end
    copyto!(scratch.ABA_patch_mat.nzval, scratch.nzval_scratch)

    numeric_refactor!(scratch.cache, scratch.ABA_patch_mat)

    copyto!(scratch.rhs_buf, view(injections_full, valid_ix, :))
    solve!(scratch.cache, scratch.rhs_buf)

    fill!(scratch.θ_full_buf, 0.0)
    scratch.θ_full_buf[valid_ix, :] .= scratch.rhs_buf

    _reconstruct_flows!(
        scratch.flow_buf,
        base_data.aux_network_matrix,
        scratch.θ_full_buf,
        d.arc_indices,
        d.delta_b,
        fb_full,
        tb_full,
        arc_susceptances,
    )

    _store_contingency_slab!(
        result, k, base_data, scratch.θ_full_buf, scratch.flow_buf, fb_full, tb_full,
    )
    # Islanding detection: connectivity BFS on the post-outage arc set.
    n_isl = _count_post_contingency_islands(base_data, mod, arc_susceptances)
    result.n_islands[k] = n_isl
    result.islanded[k] = n_isl > 1
    return
end

"""
Solve one contingency via ABA-domain Woodbury (Phase 2).

Takes the BASE KLU factorization (never refactored) and computes θ_k as
`θ_base - Z · W⁻¹ · (Uᵀ · θ_base)`. Base θ must be precomputed once outside
the loop and passed in via `θ_base` (same for all contingencies with the same
injections).
"""
function _solve_one_contingency_woodbury!(
    result::TimeContingencyPowerFlowData,
    k::Int,
    mod::PNM.NetworkModification,
    base_data::ABAPowerFlowData,
    aba::PNM.ABA_Matrix,
    θ_base::Matrix{Float64},        # n_valid × n_ts, read-only shared
    scratch::ContingencyScratch,
    fb_full::Vector{Int},
    tb_full::Vector{Int},
    arc_susceptances::Vector{Float64},
    valid_ix::AbstractVector{Int},
    woodbury_max_rank::Int,
)
    d = PNM.compute_aba_delta(aba, mod)
    M = size(d.U, 2)
    if M == 0
        fill!(scratch.θ_full_buf, 0.0)
        scratch.θ_full_buf[valid_ix, :] .= θ_base
        LinearAlgebra.mul!(
            scratch.flow_buf,
            transpose(base_data.aux_network_matrix.data),
            scratch.θ_full_buf,
        )
        _store_contingency_slab!(
            result, k, base_data, scratch.θ_full_buf, scratch.flow_buf,
            fb_full, tb_full,
        )
        # No arc modifications → base connectivity unchanged.
        result.n_islands[k] = 1
        result.islanded[k] = false
        return
    end
    if M > woodbury_max_rank
        throw(
            ArgumentError(
                "Woodbury requested but M=$M exceeds woodbury_max_rank=$woodbury_max_rank.",
            ),
        )
    end

    n_valid = size(aba.data, 1)
    Z_buf = Matrix{Float64}(undef, n_valid, M)
    # `scratch.base_K` is per-task; `aba.K` would race on `klu_common` writes.
    wf = PNM.compute_aba_woodbury_factors(
        scratch.base_K, d.U, d.delta_b, d.arc_indices, Z_buf,
    )
    θ_k_reduced = copy(θ_base)
    scratch_MxT = Matrix{Float64}(undef, M, size(θ_base, 2))
    PNM.apply_aba_woodbury_correction!(θ_k_reduced, d.U, wf, scratch_MxT)

    fill!(scratch.θ_full_buf, 0.0)
    scratch.θ_full_buf[valid_ix, :] .= θ_k_reduced

    _reconstruct_flows!(
        scratch.flow_buf,
        base_data.aux_network_matrix,
        scratch.θ_full_buf,
        d.arc_indices,
        d.delta_b,
        fb_full,
        tb_full,
        arc_susceptances,
    )

    _store_contingency_slab!(
        result, k, base_data, scratch.θ_full_buf, scratch.flow_buf,
        fb_full, tb_full,
    )
    # Islanding report: graph BFS is the ground truth. Cross-check against
    # Woodbury's numerical detector: if Woodbury flagged islanding but the
    # graph says connected, the W matrix was just ill-conditioned (not a true
    # topology split) — trust the graph.
    n_isl = _count_post_contingency_islands(base_data, mod, arc_susceptances)
    result.n_islands[k] = n_isl
    result.islanded[k] = n_isl > 1
    return
end

"""Stamp θ / flows / loss / angle-diff into 3D slabs. Does NOT touch
`result.converged` — that's a BitMatrix with 1-bit-per-cell packing that races
on the underlying word across threads. Convergence is tracked via a per-task
Vector{Bool} and merged serially after the parallel loop."""
function _store_contingency_slab!(
    result::TimeContingencyPowerFlowData,
    k::Int,
    base_data::ABAPowerFlowData,
    θ_full::Matrix{Float64},
    flow::Matrix{Float64},
    fb_full::Vector{Int},
    tb_full::Vector{Int},
)
    # Invariant fields (bus_magnitude, bus_type, injections, withdrawals)
    # are populated once across all slabs by `_populate_base_slab!` — only
    # contingency-varying fields are written here.
    @inbounds begin
        result.bus_angles[:, :, k] .= θ_full
        result.arc_active_power_flow_from_to[:, :, k] .= flow
        result.arc_active_power_flow_to_from[:, :, k] .= .-flow
        if result.arc_active_power_losses !== nothing
            Rs = _get_arc_resistances(base_data)
            result.arc_active_power_losses[:, :, k] .= Rs .* flow .^ 2
        end
    end
    _compute_arc_angle_differences_slab!(
        view(result.arc_angle_differences, :, :, k), θ_full, fb_full, tb_full,
    )
    return
end

# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

"""
    solve_contingency_dc_power_flow!(result, base_data, sys, contingencies; kwargs...)

Run DC contingency analysis in-place into a pre-built `TimeContingencyPowerFlowData`.

# Arguments
- `result::TimeContingencyPowerFlowData` — pre-allocated result container.
  Must have at least `length(contingencies) + 1` contingency slots (slot 1 is
  the base case).
- `base_data::ABAPowerFlowData` — a `PowerFlowData` configured for `DCPowerFlow`.
  `solve_power_flow!(base_data)` must have been called already (base case
  solved).
- `sys::PSY.System` — the system that `base_data` was built from.
- `contingencies::AbstractVector` — a vector of `PSY.Outage` or
  `PNM.NetworkModification`. Each element is one contingency to evaluate.

# Keyword arguments
- `strategy::Symbol = :refactor` — `:refactor` is the only supported path in
  Phase 1. `:woodbury` and `:auto` will error until Phase 2 lands.
- `on_failure::Symbol = :warn` — `:warn` catches numerical failures per
  contingency and marks `converged[:, k] = false`; `:error` rethrows.
"""
function solve_contingency_dc_power_flow!(
    result::TimeContingencyPowerFlowData,
    base_data::ABAPowerFlowData,
    sys::PSY.System,
    contingencies::AbstractVector;
    strategy::Symbol = :auto,
    woodbury_max_rank::Int = 8,
    on_failure::Symbol = :warn,
    parallel::Bool = false,
)
    if !(strategy in (:refactor, :woodbury, :auto, :auto_calibrate))
        throw(
            ArgumentError(
                "strategy must be :refactor, :woodbury, :auto, or :auto_calibrate.",
            ),
        )
    end
    n_ctg = get_n_contingencies(result)
    if length(contingencies) + 1 != n_ctg
        throw(
            ArgumentError(
                "length(contingencies)+1 ($(length(contingencies)+1)) must equal " *
                "get_n_contingencies(result) ($n_ctg); slot 1 is the base case.",
            ),
        )
    end

    aba = base_data.power_network_matrix
    ba = base_data.aux_network_matrix

    mods, labels = _resolve_modifications(ba, sys, contingencies)
    for (i, m) in enumerate(mods)
        result.network_modifications[i + 1] = m
    end

    _populate_base_slab!(result, base_data)

    # --- Shared read-only precomputation ---
    fb_full, tb_full = _arc_endpoint_indices(base_data)
    arc_sus = PNM._get_arc_susceptances(ba)
    n_ts = get_time_steps(base_data)
    n_valid = size(aba.data, 1)
    n_arcs = length(get_arc_axis(base_data))
    n_bus_full = size(result.bus_angles, 1)

    nzval_base = copy(aba.data.nzval)

    injections_full =
        base_data.bus_active_power_injections .- base_data.bus_active_power_withdrawals
    injections_full .+= base_data.bus_hvdc_net_power
    valid_ix = _valid_ix_indices(base_data)

    # Base θ computed once (reduced basis) — read-only shared across tasks
    θ_base = Matrix{Float64}(undef, n_valid, n_ts)
    copyto!(θ_base, view(injections_full, valid_ix, :))
    LinearAlgebra.ldiv!(aba.K, θ_base)

    # --- Scratch pool (one ContingencyScratch per worker) ---
    n_workers = parallel ? Threads.nthreads() : 1
    pool = Channel{ContingencyScratch}(n_workers)
    for _ in 1:n_workers
        put!(pool, ContingencyScratch(aba, n_bus_full, n_arcs, n_ts))
    end

    # :auto_calibrate times one representative contingency both ways and
    # sets the Woodbury/refactor threshold empirically (clamped to [1, 32]).
    eff_strategy, eff_woodbury_max_rank =
        if strategy === :auto_calibrate && !isempty(mods)
            _calibrate_woodbury_rank!(
                mods, base_data, aba, θ_base, arc_sus, valid_ix,
                injections_full, nzval_base, fb_full, tb_full,
                n_bus_full, n_arcs, n_ts, woodbury_max_rank,
            )
        else
            (strategy, woodbury_max_rank)
        end

    # Per-contingency convergence tracked via a Vector{Bool} (not the shared
    # BitMatrix) to avoid word-level bit-packing races under @threads.
    ctg_converged = fill(false, length(mods))

    # `parallel=false` sets `n_workers=1` upstream, making `@threads :dynamic`
    # functionally serial (one task, one cache in the pool). No separate branch.
    Threads.@threads :dynamic for k0 in 1:length(mods)
        _dispatch_one!(
            result, k0, mods, base_data, aba, ba, θ_base,
            pool, nzval_base, fb_full, tb_full, arc_sus,
            injections_full, valid_ix, eff_strategy, eff_woodbury_max_rank,
            labels, on_failure, ctg_converged,
        )
    end

    # Serial merge of per-contingency convergence into the shared BitMatrix.
    for (k0, conv) in enumerate(ctg_converged)
        result.converged[:, k0 + 1] .= conv
    end
    return result
end

"""
    _calibrate_woodbury_rank!(mods, base_data, aba, ...) -> (Symbol, Int)

Calibration path for `strategy = :auto_calibrate`. Times a single
representative contingency (the median-M one in the input set) through both
the refactor and Woodbury paths, and picks a Woodbury rank threshold from
whichever is faster per unit of M. Returns `(:auto, calibrated_rank)`.

The calibration runs against throwaway result/scratch state — it does NOT
contaminate the caller's `result` or scratch pool. The returned rank is
clamped to `[1, 32]` to avoid pathological values (e.g., if the system is
tiny and both paths take < 1 μs).

Falls back to `(:auto, requested_woodbury_max_rank)` if calibration throws
for any reason (tapped branches triggering B1's ArgumentError, empty mod
set, etc.).
"""
function _calibrate_woodbury_rank!(
    mods::Vector{PNM.NetworkModification},
    base_data::ABAPowerFlowData,
    aba::PNM.ABA_Matrix,
    θ_base::Matrix{Float64},
    arc_sus::Vector{Float64},
    valid_ix::AbstractVector{Int},
    injections_full::Matrix{Float64},
    nzval_base::Vector{Float64},
    fb_full::Vector{Int},
    tb_full::Vector{Int},
    n_bus_full::Int,
    n_arcs::Int,
    n_ts::Int,
    fallback_rank::Int,
)::Tuple{Symbol, Int}
    try
        # Pick the contingency with median M so timing reflects typical cost.
        Ms = [length(m.arc_modifications) for m in mods]
        if all(iszero, Ms)
            return (:auto, fallback_rank)
        end
        ord = sortperm(Ms)
        calib_idx = ord[clamp(div(length(mods) + 1, 2), 1, length(mods))]
        mod = mods[calib_idx]
        M = length(mod.arc_modifications)
        if M == 0
            return (:auto, fallback_rank)
        end

        tmp_result = TimeContingencyPowerFlowData(
            n_bus_full, n_arcs, n_ts, ["base", "calibration"];
            make_arc_active_power_losses = true,
        )
        _populate_base_slab!(tmp_result, base_data)
        tmp_scratch = ContingencyScratch(aba, n_bus_full, n_arcs, n_ts)

        # Warm-up pass to exclude JIT compilation from timing.
        _solve_one_contingency_refactor!(
            tmp_result, 2, mod, base_data, aba, tmp_scratch, nzval_base,
            fb_full, tb_full, arc_sus, injections_full, valid_ix,
        )
        _solve_one_contingency_woodbury!(
            tmp_result, 2, mod, base_data, aba, θ_base, tmp_scratch,
            fb_full, tb_full, arc_sus, valid_ix, max(M, fallback_rank),
        )

        # Median-of-N timing to suppress single-run jitter. N=5 trades small
        # extra startup cost (≈ 5× one solve) for stability.
        N = 5
        t_refactors = Vector{Float64}(undef, N)
        t_woodburys = Vector{Float64}(undef, N)
        for i in 1:N
            t_refactors[i] = @elapsed _solve_one_contingency_refactor!(
                tmp_result, 2, mod, base_data, aba, tmp_scratch, nzval_base,
                fb_full, tb_full, arc_sus, injections_full, valid_ix,
            )
            t_woodburys[i] = @elapsed _solve_one_contingency_woodbury!(
                tmp_result, 2, mod, base_data, aba, θ_base, tmp_scratch,
                fb_full, tb_full, arc_sus, valid_ix, max(M, fallback_rank),
            )
        end
        sort!(t_refactors)
        sort!(t_woodburys)
        t_refactor = t_refactors[(N + 1) ÷ 2]
        t_woodbury = t_woodburys[(N + 1) ÷ 2]

        # Below ~100 ns the @elapsed resolution is too coarse to draw a
        # reliable conclusion. Fall back to the requested rank in that case.
        if t_refactor < 1e-7 || t_woodbury < 1e-7
            @info "Contingency DC PF :auto_calibrate: timings below resolution; using fallback rank." fallback_rank =
                fallback_rank
            return (:auto, fallback_rank)
        end

        woodbury_per_M = t_woodbury / max(M, 1)
        estimated_crossover = floor(Int, t_refactor / woodbury_per_M)
        calibrated = clamp(estimated_crossover, 1, 32)
        @info "Contingency DC PF strategy calibration" M = M t_refactor =
            t_refactor t_woodbury = t_woodbury calibrated_woodbury_max_rank =
            calibrated
        return (:auto, calibrated)
    catch e
        if e isa ArgumentError || e isa LinearAlgebra.SingularException ||
           e isa DomainError
            @warn "Contingency DC PF :auto_calibrate failed numerically; using requested woodbury_max_rank." exception =
                (e, catch_backtrace()) fallback_rank = fallback_rank
            return (:auto, fallback_rank)
        else
            rethrow()
        end
    end
end

"""Dispatch one contingency: take scratch from pool, pick strategy, solve,
return scratch. Failure isolation wrapped here."""
function _dispatch_one!(
    result, k0, mods, base_data, aba, ba, θ_base,
    pool, nzval_base, fb_full, tb_full, arc_sus,
    injections_full, valid_ix, strategy, woodbury_max_rank,
    labels, on_failure, ctg_converged,
)
    k = k0 + 1
    mod = mods[k0]
    M = length(mod.arc_modifications)
    eff_strategy = strategy
    if strategy === :auto
        eff_strategy = (M <= woodbury_max_rank) ? :woodbury : :refactor
    elseif strategy === :woodbury && M > woodbury_max_rank
        @warn "Contingency $(labels[k0]): M=$M > woodbury_max_rank=$woodbury_max_rank; falling back to :refactor."
        eff_strategy = :refactor
    end

    scratch = take!(pool)
    try
        try
            if eff_strategy === :woodbury
                _solve_one_contingency_woodbury!(
                    result, k, mod, base_data, aba, θ_base, scratch,
                    fb_full, tb_full, arc_sus, valid_ix, woodbury_max_rank,
                )
            else
                _solve_one_contingency_refactor!(
                    result, k, mod, base_data, aba, scratch, nzval_base,
                    fb_full, tb_full, arc_sus, injections_full, valid_ix,
                )
            end
            # Thread-safe write: each index is a distinct byte in Vector{Bool}.
            ctg_converged[k0] = true
        catch e
            if _is_numerical_failure(e) && on_failure === :warn
                @warn "Contingency $(labels[k0]) failed numerically; marking not converged." exception =
                    (e, catch_backtrace())
                ctg_converged[k0] = false
            else
                rethrow()
            end
        end
    finally
        put!(pool, scratch)
    end
    return
end

"""
    solve_contingency_dc_power_flow(sys, contingencies; time_steps, kwargs...)

Convenience entry point: builds `ABAPowerFlowData`, solves the base case, then
runs the contingency loop. Returns a populated `TimeContingencyPowerFlowData`.
"""
function solve_contingency_dc_power_flow(
    sys::PSY.System,
    contingencies::AbstractVector;
    time_steps::Int = 1,
    time_step_names::Vector{String} = String[],
    network_reductions::Vector{PNM.NetworkReduction} = PNM.NetworkReduction[],
    strategy::Symbol = :auto,
    woodbury_max_rank::Int = 8,
    on_failure::Symbol = :warn,
    parallel::Bool = false,
)
    pf = DCPowerFlow(;
        time_steps = time_steps,
        time_step_names = time_step_names,
        network_reductions = network_reductions,
    )
    base_data = PowerFlowData(pf, sys)
    solve_power_flow!(base_data)

    n_buses = size(base_data.bus_angles, 1)
    n_arcs = size(base_data.arc_active_power_flow_from_to, 1)
    n_lccs = get_lcc_count(base_data)
    n_ctg = length(contingencies) + 1
    labels = String["base"]
    # Stamp user-facing labels now (the internal loop repeats this; cheap)
    for (i, c) in enumerate(contingencies)
        lab = c isa PSY.Outage ? _contingency_label(c, i) :
              (isempty(c.label) ? "ctg_$(i)" : c.label)
        push!(labels, lab)
    end
    result = TimeContingencyPowerFlowData(
        n_buses,
        n_arcs,
        time_steps,
        labels;
        n_lccs = n_lccs,
        make_arc_active_power_losses = true,
    )
    solve_contingency_dc_power_flow!(
        result,
        base_data,
        sys,
        contingencies;
        strategy = strategy,
        woodbury_max_rank = woodbury_max_rank,
        on_failure = on_failure,
        parallel = parallel,
    )
    return result
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""Return an explicit integer vector of non-ref bus positions.

`get_valid_ix(base_data)` returns a `Not{Vector{Int}}` (an `InvertedIndex`) that
works for indexing but is inconvenient for explicit row-copying from a full-axis
matrix into a reduced-axis matrix. This helper returns the complementary index
set as a concrete `Vector{Int}`."""
function _valid_ix_indices(base_data::ABAPowerFlowData)
    ref_positions = PNM.get_ref_bus_position(get_metadata_matrix(base_data))
    n = size(base_data.bus_angles, 1)
    is_ref = falses(n)
    for p in ref_positions
        is_ref[p] = true
    end
    return [i for i in 1:n if !is_ref[i]]
end

"""Numerical failures we swallow per contingency under `on_failure=:warn`.
Lexical bugs (MethodError, BoundsError, DimensionMismatch) rethrow.

KLU throws `LinearAlgebra.SingularException` on zero pivots (via `kluerror`);
numeric_refactor! throws `ArgumentError("different sparse structure")` if the
pattern drifts."""
_is_numerical_failure(e::LinearAlgebra.SingularException) = true
_is_numerical_failure(e::DomainError) = true
_is_numerical_failure(e::ArgumentError) =
    occursin("different sparse structure", sprint(showerror, e))
_is_numerical_failure(e) = false
