"""
    jacobian_indefiniteness_detection.jl

Functions for detecting indefiniteness of Jacobian matrices before Newton-Raphson solver.
Employs sparse factorization methods to efficiently determine the eigenvalue spectrum
without expensive dense computations.

Key Methods:
- Inertia computation via sparse LU factorization (exploits sparsity)
- Quick definiteness checks (positive definite, negative definite, indefinite)
- Eigenvalue bounds and spectrum classification
"""

"""
    InertiaResult

Result of inertia computation containing eigenvalue signature.

Fields:
- `n_positive::Int`: Number of positive eigenvalues
- `n_negative::Int`: Number of negative eigenvalues
- `n_zero::Int`: Number of zero eigenvalues (within tolerance)
- `is_indefinite::Bool`: True if both positive and negative eigenvalues exist
- `is_positive_definite::Bool`: True if all eigenvalues are positive
- `is_negative_definite::Bool`: True if all eigenvalues are negative
- `tolerance::Float64`: Tolerance used for zero eigenvalue detection
"""
struct InertiaResult
    n_positive::Int
    n_negative::Int
    n_zero::Int
    is_indefinite::Bool
    is_positive_definite::Bool
    is_negative_definite::Bool
    tolerance::Float64
end

function Base.show(io::IO, ir::InertiaResult)
    println(io, "InertiaResult:")
    println(io, "  Positive eigenvalues: $(ir.n_positive)")
    println(io, "  Negative eigenvalues: $(ir.n_negative)")
    println(io, "  Zero eigenvalues (tol=$(ir.tolerance)): $(ir.n_zero)")
    println(io, "  Indefinite: $(ir.is_indefinite)")
    println(io, "  Positive Definite: $(ir.is_positive_definite)")
    println(io, "  Negative Definite: $(ir.is_negative_definite)")
end

"""
    compute_inertia_via_sparse_lu(J::SparseMatrixCSC{Float64,Int32};
                                   tolerance::Float64=1e-14) -> InertiaResult

Compute the inertia (eigenvalue signature) of a sparse matrix using sparse LU factorization.

This method exploits sparsity by leveraging the KLU sparse factorization. The inertia is
estimated from the signs of pivot elements in the LU decomposition. This approach is much
faster than dense eigenvalue computation, especially for large sparse matrices.

# Arguments
- `J::SparseMatrixCSC{Float64,Int32}`: Sparse Jacobian matrix
- `tolerance::Float64`: Tolerance for determining zero eigenvalues (default: 1e-14)

# Returns
- `InertiaResult`: Structure containing inertia information and definiteness status

# Theoretical Background
The Sylvester's law of inertia states that the inertia of a matrix is preserved under
congruence transformations. The LU factorization provides an approximate inertia by:
1. Counting positive pivots (relate to positive eigenvalues)
2. Counting negative pivots (relate to negative eigenvalues)
3. Detecting small/zero pivots (relate to near-zero eigenvalues)

Note: This is an approximation for non-symmetric matrices. For the power flow Jacobian
(non-symmetric), we examine both the pure factorization and the symmetric part.

# Example
```julia
inertia = compute_inertia_via_sparse_lu(J)
if inertia.is_indefinite
    @warn "Jacobian is indefinite - Newton-Raphson may have convergence issues"
end
```
"""
function compute_inertia_via_sparse_lu(J::SparseMatrixCSC{Float64,Int32};
                                       tolerance::Float64=1e-14)
    n = size(J, 1)

    try
        # Attempt LU factorization using KLU for sparse matrices
        F = KLU.klu(J)

        # Extract diagonal elements from U (the upper triangular factor)
        # For LU decomposition J = P*L*U*Q, the diagonal of U contains the "pivots"
        diag_U = diag(F.U)

        # Count positive and negative pivots
        n_positive = count(x -> x > tolerance, diag_U)
        n_negative = count(x -> x < -tolerance, diag_U)
        n_zero = count(x -> abs(x) <= tolerance, diag_U)

        is_indefinite = (n_positive > 0) && (n_negative > 0)
        is_positive_definite = (n_negative == 0) && (n_zero == 0)
        is_negative_definite = (n_positive == 0) && (n_zero == 0)

        return InertiaResult(
            n_positive, n_negative, n_zero,
            is_indefinite, is_positive_definite, is_negative_definite,
            tolerance
        )
    catch e
        # If LU factorization fails, the matrix is likely singular or severely ill-conditioned
        # This can indicate indefiniteness or numerical issues
        @warn "LU factorization failed: $(e). Matrix may be singular or ill-conditioned."
        return InertiaResult(0, 0, 0, true, false, false, tolerance)
    end
end

"""
    is_jacobian_indefinite(J::SparseMatrixCSC{Float64,Int32};
                          tolerance::Float64=1e-14) -> Bool

Quick check to determine if a Jacobian matrix is indefinite.

Returns `true` if the matrix has both positive and negative eigenvalues, indicating
indefiniteness. This is a fast check that exploits sparsity via LU factorization.

# Arguments
- `J::SparseMatrixCSC{Float64,Int32}`: Sparse Jacobian matrix
- `tolerance::Float64`: Tolerance for zero eigenvalue detection

# Returns
- `Bool`: `true` if indefinite, `false` otherwise
"""
function is_jacobian_indefinite(J::SparseMatrixCSC{Float64,Int32};
                               tolerance::Float64=1e-14)
    inertia = compute_inertia_via_sparse_lu(J; tolerance=tolerance)
    return inertia.is_indefinite
end

"""
    is_positive_definite(J::SparseMatrixCSC{Float64,Int32};
                        tolerance::Float64=1e-14) -> Bool

Check if a sparse matrix is positive definite.

A matrix is positive definite if all its eigenvalues are positive.

# Arguments
- `J::SparseMatrixCSC{Float64,Int32}`: Sparse matrix
- `tolerance::Float64`: Tolerance for zero eigenvalue detection

# Returns
- `Bool`: `true` if positive definite, `false` otherwise
"""
function is_positive_definite(J::SparseMatrixCSC{Float64,Int32};
                             tolerance::Float64=1e-14)
    inertia = compute_inertia_via_sparse_lu(J; tolerance=tolerance)
    return inertia.is_positive_definite
end

"""
    is_negative_definite(J::SparseMatrixCSC{Float64,Int32};
                        tolerance::Float64=1e-14) -> Bool

Check if a sparse matrix is negative definite.

A matrix is negative definite if all its eigenvalues are negative.

# Arguments
- `J::SparseMatrixCSC{Float64,Int32}`: Sparse matrix
- `tolerance::Float64`: Tolerance for zero eigenvalue detection

# Returns
- `Bool`: `true` if negative definite, `false` otherwise
"""
function is_negative_definite(J::SparseMatrixCSC{Float64,Int32};
                             tolerance::Float64=1e-14)
    inertia = compute_inertia_via_sparse_lu(J; tolerance=tolerance)
    return inertia.is_negative_definite
end

"""
    check_jacobian_symmetric_part(J::SparseMatrixCSC{Float64,Int32};
                                  tolerance::Float64=1e-14) -> InertiaResult

Analyze the symmetric part of a non-symmetric Jacobian matrix.

For a non-symmetric Jacobian J, this function computes the inertia of the symmetric
part: (J + J^T) / 2. This provides insight into the "symmetrized" eigenvalue spectrum,
which can be useful for understanding the local behavior near the solution.

# Arguments
- `J::SparseMatrixCSC{Float64,Int32}`: Sparse Jacobian matrix
- `tolerance::Float64`: Tolerance for zero eigenvalue detection

# Returns
- `InertiaResult`: Inertia of the symmetric part
"""
function check_jacobian_symmetric_part(J::SparseMatrixCSC{Float64,Int32};
                                       tolerance::Float64=1e-14)
    # Compute symmetric part: (J + J^T) / 2
    # Use sparse operations to maintain efficiency
    J_symmetric = (J + J') ./ 2

    return compute_inertia_via_sparse_lu(J_symmetric; tolerance=tolerance)
end

"""
    quick_indefiniteness_check(jacobian::ACPowerFlowJacobian;
                              tolerance::Float64=1e-14) -> Bool

Convenience function for quick indefiniteness check on ACPowerFlowJacobian objects.

This function provides a fast check directly on the ACPowerFlowJacobian structure
used in PowerFlows.jl, extracting the sparse matrix and checking indefiniteness.

# Arguments
- `jacobian::ACPowerFlowJacobian`: Power flow Jacobian object
- `tolerance::Float64`: Tolerance for zero eigenvalue detection

# Returns
- `Bool`: `true` if indefinite, `false` otherwise

# Example
```julia
jacobian = ACPowerFlowJacobian(data, time_step)
if quick_indefiniteness_check(jacobian)
    @warn "Jacobian is indefinite at this iteration"
    # Consider alternative solver or regularization
end
```
"""
function quick_indefiniteness_check(jacobian::ACPowerFlowJacobian;
                                    tolerance::Float64=1e-14)
    return is_jacobian_indefinite(jacobian.Jv; tolerance=tolerance)
end

"""
    get_inertia_report(J::SparseMatrixCSC{Float64,Int32};
                      tolerance::Float64=1e-14) -> String

Generate a detailed text report of matrix inertia properties.

Useful for debugging and understanding Jacobian behavior during solver iterations.

# Arguments
- `J::SparseMatrixCSC{Float64,Int32}`: Sparse matrix
- `tolerance::Float64`: Tolerance for zero eigenvalue detection

# Returns
- `String`: Formatted report of inertia properties
"""
function get_inertia_report(J::SparseMatrixCSC{Float64,Int32};
                           tolerance::Float64=1e-14)
    inertia = compute_inertia_via_sparse_lu(J; tolerance=tolerance)
    sym_inertia = check_jacobian_symmetric_part(J; tolerance=tolerance)

    n = size(J, 1)
    nnz = nnz(J)
    sparsity = 1.0 - nnz / (n * n)

    report = """
    ═══════════════════════════════════════════════════════════
    Jacobian Matrix Inertia Report
    ═══════════════════════════════════════════════════════════
    Matrix Dimensions: $(n) × $(n)
    Non-zeros: $(nnz) (sparsity: $(round(sparsity*100; digits=2))%)

    Full Jacobian J:
    ────────────────────────────────────────────────────────────
      Positive eigenvalues: $(inertia.n_positive)
      Negative eigenvalues: $(inertia.n_negative)
      Zero eigenvalues (tol=$(tolerance)): $(inertia.n_zero)
      Status: $(inertia.is_indefinite ? "INDEFINITE ⚠" :
              (inertia.is_positive_definite ? "Positive Definite ✓" :
              (inertia.is_negative_definite ? "Negative Definite ✓" : "Singular")))

    Symmetric Part (J + J^T)/2:
    ────────────────────────────────────────────────────────────
      Positive eigenvalues: $(sym_inertia.n_positive)
      Negative eigenvalues: $(sym_inertia.n_negative)
      Zero eigenvalues (tol=$(tolerance)): $(sym_inertia.n_zero)
      Status: $(sym_inertia.is_indefinite ? "INDEFINITE ⚠" :
              (sym_inertia.is_positive_definite ? "Positive Definite ✓" :
              (sym_inertia.is_negative_definite ? "Negative Definite ✓" : "Singular")))
    ═══════════════════════════════════════════════════════════
    """

    return report
end

"""
    monitor_jacobian_definiteness(jacobian::ACPowerFlowJacobian;
                                 verbose::Bool=true) -> InertiaResult

Monitor and report the definiteness properties of a Jacobian during solving.

This function provides detailed feedback useful during Newton-Raphson iterations
for diagnosing convergence issues related to Jacobian conditioning.

# Arguments
- `jacobian::ACPowerFlowJacobian`: Power flow Jacobian object
- `verbose::Bool`: If true, print detailed report

# Returns
- `InertiaResult`: Full inertia information
"""
function monitor_jacobian_definiteness(jacobian::ACPowerFlowJacobian;
                                       verbose::Bool=true)
    inertia = compute_inertia_via_sparse_lu(jacobian.Jv)

    if verbose
        println(get_inertia_report(jacobian.Jv))
    end

    return inertia
end
