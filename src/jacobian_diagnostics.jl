"""
    JacobianDiagnostics

Result of Jacobian matrix health diagnostics computed from KLU factorization.

# Fields
- `rcond::Float64`: Reciprocal condition number estimate (0 = singular, 1 = perfectly conditioned).
- `condest::Float64`: Estimated 1-norm condition number.
- `rgrowth::Float64`: Pivot growth factor — large values indicate numerical instability risk.
- `numerical_rank::Int`: Numerical rank detected by KLU.
- `n::Int`: Matrix dimension.
- `is_singular::Bool`: `true` if `numerical_rank < n`.
- `is_ill_conditioned::Bool`: `true` if `rcond < tolerance`.
- `min_pivot::Float64`: Minimum absolute pivot value from diag(U).
- `max_pivot::Float64`: Maximum absolute pivot value from diag(U).
"""
struct JacobianDiagnostics
    rcond::Float64
    condest::Float64
    rgrowth::Float64
    numerical_rank::Int
    n::Int
    is_singular::Bool
    is_ill_conditioned::Bool
    min_pivot::Float64
    max_pivot::Float64
end

function Base.show(io::IO, d::JacobianDiagnostics)
    println(io, "JacobianDiagnostics ($(d.n)×$(d.n)):")
    println(io, "  rcond:          $(d.rcond)")
    println(io, "  condest:        $(d.condest)")
    println(io, "  rgrowth:        $(d.rgrowth)")
    println(io, "  numerical_rank: $(d.numerical_rank) / $(d.n)")
    println(io, "  pivot range:    [$(d.min_pivot), $(d.max_pivot)]")
    status = d.is_singular ? "SINGULAR" : (d.is_ill_conditioned ? "ILL-CONDITIONED" : "OK")
    print(io, "  status:         $(status)")
end

"""
    _call_klu_diagnostics(K::KLU.KLUFactorization{Float64, T},
                          rf_common, A_nzval) where {T <: TIs}

Call KLU diagnostic functions (rcond, condest, rgrowth) on an existing factorization.
Returns `(rcond, condest, rgrowth)`. On failure of any call, the corresponding value is NaN.
"""
function _call_klu_diagnostics(
    K::KLU.KLUFactorization{Float64, T},
    rf_common::Union{Base.RefValue{KLU.klu_common}, Base.RefValue{KLU.klu_l_common}},
    A_nzval::Vector{Float64},
) where {T <: TIs}
    sym_ptr = if T === Int32
        Ptr{KLU.klu_symbolic}(K._symbolic)
    else
        Ptr{KLU.klu_l_symbolic}(K._symbolic)
    end
    num_ptr =
        T === Int32 ? Ptr{KLU.klu_numeric}(K._numeric) :
        Ptr{KLU.klu_l_numeric}(K._numeric)

    rcond_fn = T === Int32 ? KLU.klu_rcond : KLU.klu_l_rcond
    condest_fn = T === Int32 ? KLU.klu_condest : KLU.klu_l_condest
    rgrowth_fn = T === Int32 ? KLU.klu_rgrowth : KLU.klu_l_rgrowth

    rc = rcond_fn(sym_ptr, num_ptr, rf_common) == 1 ? K.common.rcond : NaN
    ce = if condest_fn(K.colptr, A_nzval, sym_ptr, num_ptr, rf_common) == 1
        K.common.condest
    else
        NaN
    end
    rg = if rgrowth_fn(K.colptr, K.rowval, A_nzval, sym_ptr, num_ptr, rf_common) == 1
        K.common.rgrowth
    else
        NaN
    end

    return (rc, ce, rg)
end

"""
    compute_jacobian_diagnostics(cache::KLULinSolveCache{T},
                                 A::SparseMatrixCSC{Float64, T};
                                 tolerance::Float64=DEFAULT_RCOND_TOLERANCE
                                 ) where {T <: TIs} -> JacobianDiagnostics

Compute diagnostics for a Jacobian matrix using an existing KLU factorization cache.

This is the preferred method during solver iterations — the factorization already exists
in the cache, so the only extra cost is the O(n) diagnostic calls.

`A` must be the same matrix that was most recently factored into `cache`.
"""
function compute_jacobian_diagnostics(
    cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T};
    tolerance::Float64 = DEFAULT_RCOND_TOLERANCE,
) where {T <: TIs}
    K = cache.K
    n = K.n
    rcond_val, condest_val, rgrowth_val =
        _call_klu_diagnostics(K, cache.rf_common, A.nzval)
    numerical_rank = K.common.numerical_rank

    U_diag = LinearAlgebra.diag(K.U)
    abs_pivots = abs.(U_diag)
    min_piv = minimum(abs_pivots)
    max_piv = maximum(abs_pivots)

    is_sing = numerical_rank < n
    is_ill = !isnan(rcond_val) && rcond_val < tolerance

    return JacobianDiagnostics(
        rcond_val, condest_val, rgrowth_val,
        numerical_rank, n,
        is_sing, is_ill,
        min_piv, max_piv,
    )
end

"""
    compute_jacobian_diagnostics(J::SparseMatrixCSC{Float64, T};
                                 tolerance::Float64=DEFAULT_RCOND_TOLERANCE
                                 ) where {T <: TIs} -> JacobianDiagnostics

Compute diagnostics for a sparse matrix by creating a temporary KLU factorization.

Use this for one-off checks outside the solver loop (e.g., pre-solve validation).
"""
function compute_jacobian_diagnostics(
    J::SparseMatrixCSC{Float64, T};
    tolerance::Float64 = DEFAULT_RCOND_TOLERANCE,
) where {T <: TIs}
    n = size(J, 1)
    cache = KLULinSolveCache(J)
    try
        full_factor!(cache, J)
    catch e
        if e isa LinearAlgebra.SingularException
            return JacobianDiagnostics(
                0.0, Inf, Inf,
                cache.K.common.numerical_rank, n,
                true, true,
                0.0, 0.0,
            )
        end
        rethrow(e)
    end
    return compute_jacobian_diagnostics(cache, J; tolerance = tolerance)
end

"""
    quick_jacobian_check(cache::KLULinSolveCache{T};
                         tolerance::Float64=DEFAULT_RCOND_TOLERANCE
                         ) where {T <: TIs} -> Bool

Lightweight health check using only rcond. Returns `true` if the Jacobian is non-singular
and reasonably well-conditioned (`rcond >= tolerance`).
"""
function quick_jacobian_check(
    cache::KLULinSolveCache{T};
    tolerance::Float64 = DEFAULT_RCOND_TOLERANCE,
) where {T <: TIs}
    K = cache.K
    sym_ptr = if T === Int32
        Ptr{KLU.klu_symbolic}(K._symbolic)
    else
        Ptr{KLU.klu_l_symbolic}(K._symbolic)
    end
    num_ptr =
        T === Int32 ? Ptr{KLU.klu_numeric}(K._numeric) :
        Ptr{KLU.klu_l_numeric}(K._numeric)

    rcond_fn = T === Int32 ? KLU.klu_rcond : KLU.klu_l_rcond
    ok = rcond_fn(sym_ptr, num_ptr, cache.rf_common)
    if ok != 1
        return false
    end
    return K.common.rcond >= tolerance
end

"""
    quick_jacobian_check(jacobian::ACPowerFlowJacobian;
                         tolerance::Float64=DEFAULT_RCOND_TOLERANCE) -> Bool

Convenience wrapper that creates a temporary factorization and checks health.
"""
function quick_jacobian_check(
    jacobian::ACPowerFlowJacobian;
    tolerance::Float64 = DEFAULT_RCOND_TOLERANCE,
)
    cache = KLULinSolveCache(jacobian.Jv)
    full_factor!(cache, jacobian.Jv)
    return quick_jacobian_check(cache; tolerance = tolerance)
end

"""
    get_diagnostics_report(diag::JacobianDiagnostics) -> String

Generate a formatted text report of Jacobian diagnostics for debugging.
"""
function get_diagnostics_report(diag::JacobianDiagnostics)
    status =
        diag.is_singular ? "SINGULAR" :
        (diag.is_ill_conditioned ? "ILL-CONDITIONED" : "OK")
    return """
    Jacobian Diagnostics Report ($(diag.n)x$(diag.n))
    ────────────────────────────────────────────
      Reciprocal condition number: $(diag.rcond)
      Condition number estimate:   $(diag.condest)
      Pivot growth factor:         $(diag.rgrowth)
      Numerical rank:              $(diag.numerical_rank) / $(diag.n)
      Pivot range:                 [$(diag.min_pivot), $(diag.max_pivot)]
      Status:                      $(status)
    ────────────────────────────────────────────"""
end
