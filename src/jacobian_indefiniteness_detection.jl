"""
Voltage stability diagnostics via reduced Jacobian eigenvalues.

The reduced Jacobian J_R = J_QV - J_Qθ · J_Pθ⁻¹ · J_PV relates to voltage stability:
positive definiteness of J_R means the system is voltage-stable. Loss of positive
definiteness signals a voltage collapse bifurcation.

Uses ARPACK (implicitly restarted Lanczos) to compute extreme eigenvalues of J_R
without materializing the dense reduced Jacobian. The operator v -> J_R v is applied
via sparse sub-block multiplies and a KLU factorization of J_Pθ.
"""

"""
    DefinitenessResult

Result of eigenvalue-based definiteness analysis on the reduced Jacobian.
"""
struct DefinitenessResult
    is_positive_definite::Bool
    is_negative_definite::Bool
    smallest_eigenvalue::Float64
    largest_eigenvalue::Float64
    nearest_zero_eigenvalue::Float64
end

is_indefinite(r::DefinitenessResult) = !r.is_positive_definite && !r.is_negative_definite

function Base.show(io::IO, r::DefinitenessResult)
    status = if r.is_positive_definite
        "POSITIVE DEFINITE"
    elseif r.is_negative_definite
        "NEGATIVE DEFINITE"
    else
        "INDEFINITE"
    end
    print(
        io,
        "DefinitenessResult: $(status) (λ_min=$(r.smallest_eigenvalue), λ_max=$(r.largest_eigenvalue), λ_nearest_zero=$(r.nearest_zero_eigenvalue))",
    )
end

"""
    ReducedJacobianOperator

Matrix-free linear operator representing the symmetric part of the reduced Jacobian:
    S = (J_R + J_R') / 2,  where  J_R = J_QV - J_Qθ · J_Pθ⁻¹ · J_PV.

Applies `v -> S v = (J_R v + J_R' v) / 2` using sparse sub-block multiplies and
a single KLU factorization of J_Pθ (transposed solves via `transpose(F)`).
"""
struct ReducedJacobianOperator
    J_QV::SparseMatrixCSC{Float64, J_INDEX_TYPE}
    J_Qθ::SparseMatrixCSC{Float64, J_INDEX_TYPE}
    J_PV::SparseMatrixCSC{Float64, J_INDEX_TYPE}
    # Transposed sub-blocks stored explicitly for fast CSC column access
    J_QV_t::SparseMatrixCSC{Float64, J_INDEX_TYPE}
    J_Qθ_t::SparseMatrixCSC{Float64, J_INDEX_TYPE}
    J_PV_t::SparseMatrixCSC{Float64, J_INDEX_TYPE}
    F_Pθ::KLU.KLUFactorization{Float64, J_INDEX_TYPE}
    n::Int
    n_pvpq::Int
    # work buffers
    tmp_pvpq::Vector{Float64}   # length n_pvpq
    tmp_pvpq2::Vector{Float64}  # length n_pvpq, for transpose path
    tmp_pq::Vector{Float64}     # length n_pq
end

Base.size(op::ReducedJacobianOperator) = (op.n, op.n)
Base.size(op::ReducedJacobianOperator, d::Integer) = d <= 2 ? op.n : 1
Base.eltype(::ReducedJacobianOperator) = Float64
LinearAlgebra.issymmetric(::ReducedJacobianOperator) = true
LinearAlgebra.checksquare(op::ReducedJacobianOperator) = op.n

function LinearAlgebra.mul!(
    y::AbstractVector{Float64},
    op::ReducedJacobianOperator,
    x::AbstractVector{Float64},
)
    # Compute J_R v:  y = J_QV * x - J_Qθ * (J_Pθ \ (J_PV * x))
    mul!(op.tmp_pvpq, op.J_PV, x)
    ldiv!(op.F_Pθ, op.tmp_pvpq)
    mul!(y, op.J_QV, x)
    mul!(y, op.J_Qθ, op.tmp_pvpq, -1.0, 1.0)

    # Compute J_R' v:  tmp_pq = J_QV' * x - J_PV' * (J_Pθ' \ (J_Qθ' * x))
    mul!(op.tmp_pvpq2, op.J_Qθ_t, x)
    ldiv!(transpose(op.F_Pθ), op.tmp_pvpq2)
    mul!(op.tmp_pq, op.J_QV_t, x)
    mul!(op.tmp_pq, op.J_PV_t, op.tmp_pvpq2, -1.0, 1.0)

    # y = (J_R v + J_R' v) / 2
    y .= (y .+ op.tmp_pq) ./ 2
    return y
end

"""
    ReducedJacobianCache

Pre-allocated workspace for extracting sub-blocks of the interleaved full Jacobian
and building the [`ReducedJacobianOperator`](@ref). Reuse across iterations.
"""
mutable struct ReducedJacobianCache
    P_row_mask::BitVector
    Q_row_mask::BitVector
    θ_col_mask::BitVector
    V_col_mask::BitVector
    n_pvpq::Int
    n_pq::Int
end

function ReducedJacobianCache(data::ACPowerFlowData, time_step::Integer)
    bus_types = @view data.bus_type[:, time_step]
    n_bus = size(data.bus_type, 1)
    pvpq_mask = bus_types .!= (PSY.ACBusTypes.REF,)
    pq_mask = bus_types .== (PSY.ACBusTypes.PQ,)

    n_pvpq = count(pvpq_mask)
    n_pq = count(pq_mask)

    odd = repeat([true, false]; outer = n_bus)
    even = repeat([false, true]; outer = n_bus)
    pvpq_interleaved = repeat(pvpq_mask; inner = 2)
    pq_interleaved = repeat(pq_mask; inner = 2)

    return ReducedJacobianCache(
        pvpq_interleaved .& odd,    # P rows (pvpq)
        pq_interleaved .& even,     # Q rows (pq)
        pvpq_interleaved .& even,   # θ cols (pvpq)
        pq_interleaved .& odd,      # V cols (pq)
        n_pvpq,
        n_pq,
    )
end

"""
    build_reduced_jacobian_operator(cache, Jv) -> ReducedJacobianOperator

Build the matrix-free operator for J_R from the full interleaved Jacobian.
"""
function build_reduced_jacobian_operator(
    cache::ReducedJacobianCache,
    Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE},
)
    J_Pθ = Jv[cache.P_row_mask, cache.θ_col_mask]
    J_PV = Jv[cache.P_row_mask, cache.V_col_mask]
    J_Qθ = Jv[cache.Q_row_mask, cache.θ_col_mask]
    J_QV = Jv[cache.Q_row_mask, cache.V_col_mask]

    F_Pθ = KLU.klu(J_Pθ)

    return ReducedJacobianOperator(
        J_QV, J_Qθ, J_PV,
        SparseMatrixCSC{Float64, J_INDEX_TYPE}(J_QV'),
        SparseMatrixCSC{Float64, J_INDEX_TYPE}(J_Qθ'),
        SparseMatrixCSC{Float64, J_INDEX_TYPE}(J_PV'),
        F_Pθ,
        cache.n_pq,
        cache.n_pvpq,
        Vector{Float64}(undef, cache.n_pvpq),
        Vector{Float64}(undef, cache.n_pvpq),
        Vector{Float64}(undef, cache.n_pq),
    )
end

"""
    check_definiteness(cache, Jv; nev, tol) -> DefinitenessResult

Check whether the reduced Jacobian is positive definite, negative definite,
or indefinite by computing extreme eigenvalues with ARPACK.

Returns both the classification and the actual smallest/largest eigenvalues.
"""
function check_definiteness(
    cache::ReducedJacobianCache,
    Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE};
    nev::Int = 1,
    ncv::Int = max(ceil(Int, sqrt(cache.n_pq)), 2 * nev + 1),
    ncv_sm::Int = max(min(cache.n_pq - 1, 2 * ceil(Int, sqrt(cache.n_pq))), 2 * nev + 1),
    maxiter::Int = 2000,
    maxiter_sm::Int = 3000,
    tol::Float64 = 1e-4,
)
    if cache.n_pq == 0
        return DefinitenessResult(true, false, Inf, -Inf, Inf)
    end
    if cache.n_pq <= ncv_sm + 1
        return _check_definiteness_dense(cache, Jv)
    end

    # TODO: could accept an existing ReducedJacobianOperator and reuse it
    # between calls (e.g. across Newton iterations), re-factorizing J_Pθ in place.
    op = build_reduced_jacobian_operator(cache, Jv)

    v0 = zeros(Float64, 0)

    # smallest algebraic eigenvalue — give it the same budget as :SM
    λ_min = _try_eigs(op, :SR; nev, ncv = ncv_sm, tol, maxiter = maxiter_sm, v0) do λs
        real(λs[end])
    end

    # largest algebraic eigenvalue — easy for Lanczos
    λ_max = _try_eigs(op, :LR; nev, ncv, tol, maxiter, v0) do λs
        real(λs[1])
    end

    # nearest-to-zero eigenvalue
    λ_nearest_zero = _try_eigs(op, :SM;
        nev, ncv = ncv_sm, tol, maxiter = maxiter_sm, v0) do λs
        real(λs[1])
    end

    is_pd = !isnan(λ_min) && λ_min > 0
    is_nd = !isnan(λ_max) && λ_max < 0
    return DefinitenessResult(is_pd, is_nd, λ_min, λ_max, λ_nearest_zero)
end

"""Try an ARPACK eigs call; on convergence failure return NaN instead of throwing."""
function _try_eigs(extract::Function, op, which::Symbol;
    nev, ncv, tol, maxiter, v0)
    try
        λs, = Arpack.eigs(op; nev, ncv, which, tol, maxiter, v0)
        return extract(λs)
    catch e
        if e isa Arpack.XYAUPD_Exception
            @warn "ARPACK failed to converge for which=:$which (maxiter=$maxiter, ncv=$ncv)"
            return NaN
        end
        rethrow(e)
    end
end

function _check_definiteness_dense(
    cache::ReducedJacobianCache,
    Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE},
)
    J_Pθ = Jv[cache.P_row_mask, cache.θ_col_mask]
    J_PV = Jv[cache.P_row_mask, cache.V_col_mask]
    J_Qθ = Jv[cache.Q_row_mask, cache.θ_col_mask]
    J_QV = Jv[cache.Q_row_mask, cache.V_col_mask]
    J_R = Matrix(J_QV) - Matrix(J_Qθ) * (Matrix(J_Pθ) \ Matrix(J_PV))
    S = LinearAlgebra.Symmetric((J_R + J_R') ./ 2)
    λs = LinearAlgebra.eigvals(S)
    λ_min = first(λs)
    λ_max = last(λs)
    _, idx = findmin(abs, λs)
    λ_nearest_zero = λs[idx]
    return DefinitenessResult(λ_min > 0, λ_max < 0, λ_min, λ_max, λ_nearest_zero)
end

function check_definiteness(jacobian::ACPowerFlowJacobian, time_step::Integer; kwargs...)
    cache = ReducedJacobianCache(jacobian.data, time_step)
    return check_definiteness(cache, jacobian.Jv; kwargs...)
end

"""
    get_definiteness_report(jacobian::ACPowerFlowJacobian, time_step::Integer) -> String
"""
function get_definiteness_report(jacobian::ACPowerFlowJacobian, time_step::Integer)
    cache = ReducedJacobianCache(jacobian.data, time_step)
    result = check_definiteness(cache, jacobian.Jv)
    n = size(jacobian.Jv, 1)
    definiteness = if result.is_positive_definite
        "POSITIVE DEFINITE"
    elseif result.is_negative_definite
        "NEGATIVE DEFINITE"
    else
        "INDEFINITE"
    end
    return """
    Reduced Jacobian Definiteness Report
    ────────────────────────────────────────────
      Full Jacobian:     $(n)×$(n)
      Reduced J_R:       $(cache.n_pq)×$(cache.n_pq) (PQ buses, matrix-free)
      Definiteness:      $(definiteness)
      λ_min:             $(result.smallest_eigenvalue)
      λ_max:             $(result.largest_eigenvalue)
      λ_nearest_zero:    $(result.nearest_zero_eigenvalue)
    ────────────────────────────────────────────"""
end

## Fiedler / algebraic connectivity analysis of J_Pθ

"""
    FiedlerResult

Result of Fiedler (algebraic connectivity) analysis on J_Pθ.

# Fields
- `λ2::Float64`: Fiedler value (algebraic connectivity). Small → weak connectivity → trouble.
- `λ_max::Float64`: Largest eigenvalue of J_Pθ.
- `condition_number::Float64`: `λ_max / λ2` — condition number of the reduced Laplacian.
- `line_stress::Vector{Float64}`: Per-arc stress = how much each arc's contribution to λ2
    has degraded from its nominal (flat-start) value. Sorted by `arc_ranking`.
- `arc_ranking::Vector{Int}`: Arc indices sorted by decreasing stress.
- `n_pvpq::Int`: Dimension of the reduced J_Pθ (REF buses removed).
"""
struct FiedlerResult
    λ2::Float64
    λ_max::Float64
    condition_number::Float64
    line_stress::Vector{Float64}
    arc_ranking::Vector{Int}
    n_pvpq::Int
end

function Base.show(io::IO, r::FiedlerResult)
    print(
        io,
        "FiedlerResult: λ₂=$(round(r.λ2; sigdigits=4)), " *
        "λ_max=$(round(r.λ_max; sigdigits=4)), " *
        "κ=$(round(r.condition_number; sigdigits=4)), " *
        "n=$(r.n_pvpq)",
    )
end

"""
    compute_fiedler(cache::ReducedJacobianCache, Jv, data, time_step; kwargs...) -> FiedlerResult

Compute Fiedler value and per-arc stress decomposition of J_Pθ.

Extracts J_Pθ from the full Jacobian, removes the REF bus row/column, then uses
ARPACK shift-invert mode (`sigma=0`) to find the smallest eigenvalues efficiently.
"""
function compute_fiedler(
    cache::ReducedJacobianCache,
    Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE},
    data::ACPowerFlowData,
    time_step::Integer;
    nev::Int = 2,
    ncv::Int = max(20, 2 * nev + 1),
    maxiter::Int = 1000,
    tol::Float64 = 1e-4,
)
    J_Pθ = Jv[cache.P_row_mask, cache.θ_col_mask]
    n = size(J_Pθ, 1)

    # Largest eigenvalue — standard Lanczos, easy
    λ_lr, = Arpack.eigs(J_Pθ; nev = 1, ncv = ncv, which = :LR, tol = tol,
        maxiter = maxiter, v0 = zeros(Float64, 0))
    λ_max = real(λ_lr[1])

    # Fiedler value via shift-invert around σ=0.
    # Request nev eigenvalues to skip the trivial zero eigenvalue if present.
    λs, vs = Arpack.eigs(J_Pθ; nev = nev, ncv = ncv, which = :LM,
        sigma = 0.0, tol = tol, maxiter = maxiter, v0 = zeros(Float64, 0))
    # Sort by real part; Fiedler value is the smallest nonzero one
    perm = sortperm(real.(λs))
    λ_sorted = real.(λs[perm])
    v_sorted = real.(vs[:, perm])

    # Pick the Fiedler value: smallest eigenvalue > tol (skip near-zero)
    fiedler_idx = findfirst(>(1e-8), λ_sorted)
    if isnothing(fiedler_idx)
        # All eigenvalues near zero — disconnected or singular
        return FiedlerResult(0.0, λ_max, Inf, Float64[], Int[], n)
    end
    λ2 = λ_sorted[fiedler_idx]
    v2 = v_sorted[:, fiedler_idx]
    v2 ./= norm(v2)

    # Per-arc stress decomposition
    line_stress = _compute_line_stress(J_Pθ, v2, data, cache, time_step)
    arc_ranking = sortperm(line_stress; rev = true)

    κ = λ_max / λ2
    return FiedlerResult(λ2, λ_max, κ, line_stress, arc_ranking, n)
end

"""
Compute per-arc stress: how much each arc's contribution to λ₂ has degraded
from its nominal (flat-start) value.

Uses the off-diagonal entries of J_Pθ directly as edge weights. The nominal weight
for each edge is |w_ij| at flat start (where cos(θ_ij)=1), approximated by the
susceptance from the Ybus. Since we don't have easy access to per-arc susceptance
here, we use the current |w_ij| as both the weight and the nominal baseline, and
stress = |w_ij| * (v2_i - v2_j)² measures each arc's current contribution to λ₂.
"""
function _compute_line_stress(
    J_Pθ::SparseMatrixCSC{Float64, J_INDEX_TYPE},
    v2::Vector{Float64},
    data::ACPowerFlowData,
    cache::ReducedJacobianCache,
    time_step::Integer,
)
    bus_types = @view data.bus_type[:, time_step]
    pvpq_indices = findall(!=(PSY.ACBusTypes.REF), bus_types)
    bus_lookup = get_bus_lookup(data)
    inv_bus_lookup = Dict(v => k for (k, v) in bus_lookup)

    # Map from original bus matrix index → reduced J_Pθ index
    reduced_lookup = Dict(pvpq_indices[i] => i for i in eachindex(pvpq_indices))

    arcs = get_arc_axis(data)
    n_arcs = length(arcs)
    stress = zeros(n_arcs)

    for (arc_idx, (from_bus, to_bus)) in enumerate(arcs)
        from_ix = bus_lookup[from_bus]
        to_ix = bus_lookup[to_bus]
        # Both buses must be in the reduced system (non-REF)
        ri = get(reduced_lookup, from_ix, 0)
        rj = get(reduced_lookup, to_ix, 0)
        if ri == 0 || rj == 0
            continue
        end
        # Off-diagonal of J_Pθ is the edge weight (negative for Laplacian-like)
        w_ij = J_Pθ[ri, rj]
        dv = v2[ri] - v2[rj]
        # Contribution to λ₂ (use absolute value since off-diags are typically negative)
        stress[arc_idx] = abs(w_ij) * dv^2
    end
    return stress
end

function compute_fiedler(jacobian::ACPowerFlowJacobian, time_step::Integer; kwargs...)
    cache = ReducedJacobianCache(jacobian.data, time_step)
    return compute_fiedler(cache, jacobian.Jv, jacobian.data, time_step; kwargs...)
end

"""
    get_fiedler_report(jacobian, time_step; n_top=10) -> String

Formatted report of Fiedler analysis with top stressed arcs.
"""
function get_fiedler_report(
    jacobian::ACPowerFlowJacobian,
    time_step::Integer;
    n_top::Int = 10,
)
    result = compute_fiedler(jacobian, time_step)
    arcs = get_arc_axis(jacobian.data)

    lines = String[]
    push!(lines, "Fiedler (Algebraic Connectivity) Report")
    push!(lines, "────────────────────────────────────────────")
    push!(lines, "  J_Pθ dimension:      $(result.n_pvpq)×$(result.n_pvpq)")
    push!(lines, "  Fiedler value (λ₂):  $(round(result.λ2; sigdigits=4))")
    push!(lines, "  λ_max:               $(round(result.λ_max; sigdigits=4))")
    push!(lines, "  Condition number:    $(round(result.condition_number; sigdigits=4))")
    if !isempty(result.arc_ranking)
        n_show = min(n_top, length(result.arc_ranking))
        push!(lines, "  Top $n_show stressed arcs:")
        for i in 1:n_show
            arc_idx = result.arc_ranking[i]
            arc = arcs[arc_idx]
            s = round(result.line_stress[arc_idx]; sigdigits = 4)
            push!(lines, "    $(arc[1]) → $(arc[2]): stress = $s")
        end
    end
    push!(lines, "────────────────────────────────────────────")
    return join(lines, "\n")
end

"""
    probe_negative_mode_to_csv(cache, Jv, data, time_step, csv_path; ...)

Compute the smallest algebraic eigenvalue and eigenvector of the symmetric reduced
Jacobian S = (J_R + J_R')/2 using the matrix-free `ReducedJacobianOperator`, then
dump the eigenvector to CSV with one row per PQ bus, columns: bus_number, entry_value.

Returns the eigenvalue.
"""
function probe_negative_mode_to_csv(
    cache::ReducedJacobianCache,
    Jv::SparseMatrixCSC{Float64, J_INDEX_TYPE},
    data::ACPowerFlowData,
    time_step::Integer,
    csv_path::AbstractString;
    nev::Int = 1,
    ncv::Int = max(min(cache.n_pq - 1, 2 * ceil(Int, sqrt(cache.n_pq))), 2 * nev + 1),
    maxiter::Int = 3000,
    tol::Float64 = 1e-4,
)
    op = build_reduced_jacobian_operator(cache, Jv)
    λs, vs = Arpack.eigs(op; nev = nev, ncv = ncv, which = :SR,
        tol = tol, maxiter = maxiter, v0 = zeros(Float64, 0))
    perm = sortperm(real.(λs))
    λ_min = real(λs[perm[1]])
    v = real.(vs[:, perm[1]])
    v ./= norm(v)

    # The reduced Jacobian rows/cols correspond to PQ buses (V variables).
    bus_types = @view data.bus_type[:, time_step]
    pq_indices = findall(==(PSY.ACBusTypes.PQ), bus_types)
    bus_axes = axes(data.power_network_matrix, 1)
    @assert length(pq_indices) == length(v) == cache.n_pq

    open(csv_path, "w") do io
        println(io, "row,bus_matrix_index,bus_number,variable,entry,abs_entry")
        for i in eachindex(v)
            bus_ix = pq_indices[i]
            bus_no = bus_axes[bus_ix]
            println(io, "$i,$bus_ix,$bus_no,V,$(v[i]),$(abs(v[i]))")
        end
    end
    @info "probe_negative_mode: λ_min=$λ_min  eigenvector written to $csv_path"
    return λ_min
end

"""
    monitor_jacobian_definiteness(jacobian::ACPowerFlowJacobian, time_step::Integer) -> DefinitenessResult

Run reduced Jacobian definiteness diagnostics and log the result.
Called from solver loops when `monitor_jacobian = true`.
"""
function monitor_jacobian_definiteness(jacobian::ACPowerFlowJacobian, time_step::Integer)
    result = check_definiteness(jacobian, time_step)
    @info "Reduced Jacobian definiteness: $(result)"
    return result
end
