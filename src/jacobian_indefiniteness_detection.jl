"""
    jacobian_indefiniteness_detection.jl

Heuristic pivot-sign diagnostics for Jacobian matrices. Uses sparse LU factorization
(KLU) to count positive, negative, and near-zero pivots in the U factor.

**Important:** For a general (non-symmetric) matrix, pivot signs from LU factorization
do NOT correspond to eigenvalue signs. These diagnostics are a fast heuristic for
detecting potential singularity or mixed-sign pivots, which can indicate convergence
issues in Newton-type solvers. They should not be interpreted as eigenvalue inertia
or used to claim positive/negative definiteness of non-symmetric matrices.

Key Methods:
- Pivot-sign counting via sparse LU factorization (exploits sparsity)
- Quick mixed-pivot checks (heuristic for potential convergence issues)
- Symmetric-part pivot analysis via LU on (J + J^T)/2
"""

"""
    PivotSignResult

Result of pivot-sign analysis from sparse LU factorization.

The counts refer to signs of diagonal entries in the U factor of the LU decomposition,
NOT eigenvalues. For non-symmetric matrices, pivot signs are a heuristic indicator
and do not determine definiteness.

Fields:
- `n_positive::Int`: Number of positive pivots in U
- `n_negative::Int`: Number of negative pivots in U
- `n_zero::Int`: Number of near-zero pivots in U (within tolerance)
- `tolerance::Float64`: Tolerance used for near-zero pivot detection
- `success::Bool`: True if factorization succeeded; false if it failed
"""
struct PivotSignResult
    n_positive::Int
    n_negative::Int
    n_zero::Int
    tolerance::Float64
    success::Bool
end

# Keep InertiaResult as an alias for backward compatibility
const InertiaResult = PivotSignResult

has_mixed_sign_pivots(r::PivotSignResult) = r.success && (r.n_positive > 0) && (r.n_negative > 0)
all_pivots_positive(r::PivotSignResult) = r.success && (r.n_negative == 0) && (r.n_zero == 0)
all_pivots_negative(r::PivotSignResult) = r.success && (r.n_positive == 0) && (r.n_zero == 0)

function Base.show(io::IO, r::PivotSignResult)
    println(io, "PivotSignResult (LU pivot-sign heuristic):")
    if r.success
        println(io, "  Positive pivots: $(r.n_positive)")
        println(io, "  Negative pivots: $(r.n_negative)")
        println(io, "  Near-zero pivots (tol=$(r.tolerance)): $(r.n_zero)")
        println(io, "  Mixed-sign pivots: $(has_mixed_sign_pivots(r))")
        println(io, "  All pivots positive: $(all_pivots_positive(r))")
        println(io, "  All pivots negative: $(all_pivots_negative(r))")
    else
        println(io, "  Factorization FAILED — results are unknown")
        println(io, "  Matrix is likely singular or severely ill-conditioned")
    end
end

"""
    compute_inertia_via_sparse_lu(J::SparseMatrixCSC{Float64,Int32};
                                   tolerance::Float64=1e-14) -> PivotSignResult

Count the signs of pivots (diagonal of U) from a sparse KLU factorization.

This is a fast heuristic that exploits sparsity. For symmetric matrices, pivot signs
from LU (without pivoting that breaks symmetry) relate to eigenvalue signs via
Sylvester's law of inertia. For non-symmetric matrices (like the power flow Jacobian),
pivot signs are NOT eigenvalue signs — use this only as a diagnostic indicator for
potential singularity or convergence issues.

# Arguments
- `J::SparseMatrixCSC{Float64,Int32}`: Sparse matrix
- `tolerance::Float64`: Tolerance for detecting near-zero pivots (default: 1e-14)

# Returns
- `PivotSignResult`: Pivot-sign counts with a `success` flag.
  When `success == false`, the factorization failed and counts are set to
  `n_zero == size(J,1)` (unknown pivots treated as zero).

# Example
```julia
result = compute_inertia_via_sparse_lu(J)
if !result.success
    @warn "Factorization failed — matrix may be singular"
elseif has_mixed_sign_pivots(result)
    @warn "Mixed-sign pivots detected — potential convergence issues"
end
```
"""
function compute_inertia_via_sparse_lu(J::SparseMatrixCSC{Float64, Int32};
    tolerance::Float64 = 1e-14)
    n = size(J, 1)

    try
        F = KLU.klu(J)
        U = F.U

        # Single pass over U diagonal to count pivot signs
        n_positive = 0
        n_negative = 0
        n_zero = 0
        for col in 1:n
            # Last nonzero in each column of upper triangular U is the diagonal
            col_end = U.colptr[col + 1] - 1
            if col_end >= U.colptr[col] && U.rowval[col_end] == col
                val = U.nzval[col_end]
                if val > tolerance
                    n_positive += 1
                elseif val < -tolerance
                    n_negative += 1
                else
                    n_zero += 1
                end
            else
                n_zero += 1
            end
        end

        return PivotSignResult(n_positive, n_negative, n_zero, tolerance, true)
    catch e
        @warn "LU factorization failed: $(e). Matrix may be singular or ill-conditioned."
        return PivotSignResult(0, 0, n, tolerance, false)
    end
end

"""
    is_jacobian_indefinite(J::SparseMatrixCSC{Float64,Int32};
                          tolerance::Float64=1e-14) -> Bool

Heuristic check: returns `true` if the LU pivots have mixed signs.

This does NOT determine eigenvalue indefiniteness for non-symmetric matrices.
It is a fast diagnostic for detecting potential convergence issues.
"""
function is_jacobian_indefinite(J::SparseMatrixCSC{Float64, Int32};
    tolerance::Float64 = 1e-14)
    result = compute_inertia_via_sparse_lu(J; tolerance = tolerance)
    return has_mixed_sign_pivots(result)
end

"""
    is_positive_definite(J::SparseMatrixCSC{Float64,Int32};
                        tolerance::Float64=1e-14) -> Bool

Heuristic check: returns `true` if all LU pivots are positive.

For non-symmetric matrices, all-positive pivots do NOT guarantee positive definiteness.
This is a necessary condition only for symmetric positive definite matrices with
no pivoting.
"""
function is_positive_definite(J::SparseMatrixCSC{Float64, Int32};
    tolerance::Float64 = 1e-14)
    result = compute_inertia_via_sparse_lu(J; tolerance = tolerance)
    return all_pivots_positive(result)
end

"""
    is_negative_definite(J::SparseMatrixCSC{Float64,Int32};
                        tolerance::Float64=1e-14) -> Bool

Heuristic check: returns `true` if all LU pivots are negative.

For non-symmetric matrices, all-negative pivots do NOT guarantee negative definiteness.
"""
function is_negative_definite(J::SparseMatrixCSC{Float64, Int32};
    tolerance::Float64 = 1e-14)
    result = compute_inertia_via_sparse_lu(J; tolerance = tolerance)
    return all_pivots_negative(result)
end

"""
    check_jacobian_symmetric_part(J::SparseMatrixCSC{Float64,Int32};
                                  tolerance::Float64=1e-14) -> PivotSignResult

Analyze the LU pivot signs of the symmetric part (J + J^T) / 2.

For symmetric matrices, LU pivot signs are more meaningful as a definiteness
diagnostic (though row-pivoting can still change signs). This computes
the symmetric part and runs the pivot-sign analysis on it.
"""
function check_jacobian_symmetric_part(J::SparseMatrixCSC{Float64, Int32};
    tolerance::Float64 = 1e-14)
    J_symmetric = (J + J') ./ 2
    return compute_inertia_via_sparse_lu(J_symmetric; tolerance = tolerance)
end

"""
    quick_indefiniteness_check(jacobian::ACPowerFlowJacobian;
                              tolerance::Float64=1e-14) -> Bool

Convenience wrapper: checks if the Jacobian's LU pivots have mixed signs.
"""
function quick_indefiniteness_check(jacobian::ACPowerFlowJacobian;
    tolerance::Float64 = 1e-14)
    return is_jacobian_indefinite(jacobian.Jv; tolerance = tolerance)
end

"""
    _format_pivot_sign_section(label::String, r::PivotSignResult,
                               tolerance::Float64) -> String

Format a section of the pivot-sign report for a single matrix analysis.
"""
function _format_pivot_sign_section(label::String, r::PivotSignResult,
    tolerance::Float64)
    return """
      $(label):
      ────────────────────────────────────────────────────────────
        Positive pivots: $(r.n_positive)
        Negative pivots: $(r.n_negative)
        Near-zero pivots (tol=$(tolerance)): $(r.n_zero)
        Factorization: $(r.success ? "OK" : "FAILED")
        Status: $(has_mixed_sign_pivots(r) ? "MIXED-SIGN PIVOTS" :
                (all_pivots_positive(r) ? "All pivots positive" :
                (all_pivots_negative(r) ? "All pivots negative" : "Has near-zero pivots")))
    """
end

"""
    get_inertia_report(J::SparseMatrixCSC{Float64,Int32};
                      tolerance::Float64=1e-14) -> String

Generate a text report of LU pivot-sign diagnostics for the matrix and its
symmetric part. Useful for debugging convergence issues.
"""
function get_inertia_report(J::SparseMatrixCSC{Float64, Int32};
    tolerance::Float64 = 1e-14)
    result = compute_inertia_via_sparse_lu(J; tolerance = tolerance)
    sym_result = check_jacobian_symmetric_part(J; tolerance = tolerance)
    return _format_inertia_report(J, result, sym_result, tolerance)
end

function _format_inertia_report(J::SparseMatrixCSC{Float64, Int32},
    result::PivotSignResult, sym_result::PivotSignResult,
    tolerance::Float64)
    n = size(J, 1)
    nz = SparseArrays.nnz(J)
    sparsity = 1.0 - nz / (n * n)

    return """
    ═══════════════════════════════════════════════════════════
    Jacobian Matrix Pivot-Sign Diagnostic Report
    ═══════════════════════════════════════════════════════════
    Matrix Dimensions: $(n) × $(n)
    Non-zeros: $(nz) (sparsity: $(round(sparsity*100; digits=2))%)

    $(_format_pivot_sign_section("Full Jacobian J (LU pivot signs — heuristic only)", result, tolerance))
    $(_format_pivot_sign_section("Symmetric Part (J + J^T)/2 (LU pivot signs)", sym_result, tolerance))
    ═══════════════════════════════════════════════════════════
    """
end

"""
    monitor_jacobian_definiteness(jacobian::ACPowerFlowJacobian;
                                 verbose::Bool=true) -> PivotSignResult

Run and optionally print pivot-sign diagnostics for a power flow Jacobian.
"""
function monitor_jacobian_definiteness(jacobian::ACPowerFlowJacobian;
    verbose::Bool = true)
    J = jacobian.Jv
    result = compute_inertia_via_sparse_lu(J)

    if verbose
        sym_result = check_jacobian_symmetric_part(J)
        println(_format_inertia_report(J, result, sym_result, result.tolerance))
    end

    return result
end
