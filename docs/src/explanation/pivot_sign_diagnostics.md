## Diagnosing Power Flow Convergence with Pivot-Sign Analysis

### The Problem

Newton-Raphson and related iterative solvers for the AC power flow equations depend on
the Jacobian matrix being well-conditioned at each iteration. When a system is near
voltage collapse, has islanded subnetworks, or contains modeling errors, the Jacobian
can become singular or nearly singular, causing the solver to diverge or stall.

Identifying *why* a power flow fails to converge is often harder than detecting that it
failed. A full eigenvalue decomposition of the Jacobian would reveal the trouble, but
for large-scale systems (tens of thousands of buses) this is prohibitively expensive.

### LU Pivot-Sign Heuristic

PowerFlows.jl provides a fast diagnostic based on the signs of the pivot elements
(diagonal entries of the U factor) from the sparse LU factorization already available
via KLU. The key idea:

1. Factor the Jacobian: ``PAQ = LU`` (with row and column permutations).
2. Extract ``\text{diag}(U)`` — the *pivots*.
3. Count how many pivots are positive, negative, or near-zero (within a tolerance).

This analysis runs in roughly the same time as a single Newton step (dominated by the
sparse factorization) and provides three useful signals:

| Signal | Interpretation |
|--------|---------------|
| **Near-zero pivots** | The matrix is (numerically) singular — at least one equation is redundant or a bus is electrically isolated. |
| **Mixed-sign pivots** | The matrix has both positive and negative pivots, which for symmetric matrices corresponds to indefiniteness. For the non-symmetric power flow Jacobian this is a *heuristic* warning of potential saddle-point or bifurcation behavior. |
| **All pivots same sign** | No singularity or sign-mixing detected — the Jacobian appears well-conditioned from a pivot perspective. |

### Important Caveats

For a **symmetric** matrix, the signs of LU pivots (without row exchanges that break
symmetry) are directly related to eigenvalue signs via Sylvester's law of inertia.
The power flow Jacobian is **not symmetric**, so:

- Pivot signs are **not** eigenvalue signs. A matrix can have all-positive pivots yet
  have negative eigenvalues, or vice versa.
- The labels `has_mixed_sign_pivots`, `all_pivots_positive`, and `all_pivots_negative`
  describe **pivot signs only**, not definiteness.
- For a more rigorous (but still heuristic) analysis of the symmetric part, use
  [`check_jacobian_symmetric_part`](@ref), which examines ``(J + J^T)/2``.

Despite these limitations, pivot-sign analysis is valuable in practice because:

- Near-zero pivots reliably indicate singularity regardless of symmetry.
- Mixed-sign pivots in the power flow Jacobian empirically correlate with convergence
  difficulties near voltage collapse points.
- The analysis is essentially free — it reuses the factorization the solver already
  computes.

### Usage

```julia
using PowerFlows

# After setting up a power flow and computing the Jacobian:
result = PowerFlows.compute_inertia_via_sparse_lu(jacobian.Jv)

if !result.success
    @error "Jacobian factorization failed — matrix is likely singular"
elseif result.n_zero > 0
    @warn "$(result.n_zero) near-zero pivots — check for isolated buses"
elseif result.has_mixed_sign_pivots
    @warn "Mixed-sign pivots — system may be near a bifurcation point"
end

# Print a full diagnostic report:
println(PowerFlows.get_inertia_report(jacobian.Jv))

# Quick check via the convenience wrapper:
if PowerFlows.quick_indefiniteness_check(jacobian)
    @warn "Mixed-sign pivots detected"
end
```

### Monitoring During Iteration

The [`monitor_jacobian_definiteness`](@ref) function can be called at each Newton
iteration to track how the pivot spectrum evolves:

```julia
for iter in 1:max_iter
    # ... Newton step ...
    result = PowerFlows.monitor_jacobian_definiteness(jacobian; verbose = true)
    if !result.success || result.n_zero > 0
        @warn "Jacobian becoming singular at iteration $iter"
        break
    end
end
```

### Relationship to Other Diagnostics

PowerFlows.jl also provides [`compute_jacobian_diagnostics`](@ref), which reports
condition number estimates, numerical rank, and pivot magnitude ranges. The two
diagnostic tools are complementary:

- **Jacobian diagnostics** (`JacobianDiagnostics`): focuses on *magnitude* — how
  well-conditioned is the matrix?
- **Pivot-sign diagnostics** (`PivotSignResult`): focuses on *sign structure* — are
  there mixed-sign or near-zero pivots?

Together they give a comprehensive picture of Jacobian health during power flow solving.
