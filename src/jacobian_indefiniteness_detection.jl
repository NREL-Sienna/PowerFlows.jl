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
    ncv::Int = max(40, 2 * nev + 1),
    maxiter::Int = 1000,
    tol::Float64 = 0.0,
)
    if cache.n_pq == 0
        return DefinitenessResult(true, false, Inf, -Inf, Inf)
    end
    if cache.n_pq <= ncv + 1
        return _check_definiteness_dense(cache, Jv)
    end

    # TODO: could accept an existing ReducedJacobianOperator and reuse it
    # between calls (e.g. across Newton iterations), re-factorizing J_Pθ in place.
    op = build_reduced_jacobian_operator(cache, Jv)

    v0 = zeros(Float64, 0)

    # smallest algebraic eigenvalue (most negative)
    λ_sm, = Arpack.eigs(op; nev = nev, ncv = ncv, which = :SR,
        tol = tol, maxiter = maxiter, v0 = v0)
    λ_min = real(λ_sm[end])

    # largest algebraic eigenvalue (most positive)
    λ_lg, = Arpack.eigs(op; nev = nev, ncv = ncv, which = :LR,
        tol = tol, maxiter = maxiter, v0 = v0)
    λ_max = real(λ_lg[1])

    # eigenvalue nearest to zero (smallest magnitude)
    λ_nz, = Arpack.eigs(op; nev = nev, ncv = ncv, which = :SM,
        tol = tol, maxiter = maxiter, v0 = v0)
    λ_nearest_zero = real(λ_nz[1])

    is_pd = λ_min > 0
    is_nd = λ_max < 0
    return DefinitenessResult(is_pd, is_nd, λ_min, λ_max, λ_nearest_zero)
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
