## Diagnosing Power Flow Convergence with Jacobian Diagnostics

### The Problem

Newton-Raphson and related iterative solvers for the AC power flow equations depend on
the Jacobian matrix being well-conditioned at each iteration. When a system is near
voltage collapse, has islanded subnetworks, or contains modeling errors, the Jacobian
can become singular or nearly singular, causing the solver to diverge or stall.

Identifying *why* a power flow fails to converge is often harder than detecting that it
failed. A full eigenvalue decomposition of the Jacobian would reveal the trouble, but
for large-scale systems (tens of thousands of buses) this is prohibitively expensive.

### Reduced Jacobian Definiteness

The power flow Jacobian has block structure:

```math
J = \begin{bmatrix} J_{P\theta} & J_{PV} \\ J_{Q\theta} & J_{QV} \end{bmatrix}
```

The **reduced Jacobian** is the Schur complement:

```math
J_R = J_{QV} - J_{Q\theta} \cdot J_{P\theta}^{-1} \cdot J_{PV}
```

This matrix relates reactive power to voltage magnitude for PQ buses. Its positive
definiteness has a direct physical interpretation: **positive definite J_R means the
system is voltage-stable**. Loss of positive definiteness signals proximity to a
voltage collapse bifurcation.

PowerFlows.jl checks definiteness of J_R using ARPACK (implicitly restarted Lanczos)
to compute the extreme eigenvalues of the symmetric part of J_R without materializing
the dense reduced Jacobian. The operator `v -> (J_R + J_R')v / 2` is applied via sparse
sub-block multiplies and a KLU factorization of J_Pθ. Three eigenvalues are computed:
the smallest (most negative), the largest (most positive), and the nearest to zero
(smallest magnitude), giving both a definitive classification and a quantitative
stability margin.

### KLU-based Condition Diagnostics

PowerFlows.jl also provides [`compute_jacobian_diagnostics`](@ref), which reports
condition number estimates, numerical rank, and pivot magnitude ranges from the KLU
factorization. These are complementary to the definiteness check:

- **rcond / condest**: How well-conditioned is the full Jacobian? (magnitude)
- **numerical rank**: Is the Jacobian singular or rank-deficient?
- **J_R definiteness**: Is the system voltage-stable? (sign structure)

### Usage

```julia
using PowerFlows

# Check reduced Jacobian definiteness:
result = PowerFlows.check_definiteness(jacobian, time_step)
if PowerFlows.is_indefinite(result)
    @warn "Reduced Jacobian is indefinite — voltage instability detected"
end

# Inspect eigenvalues directly:
result.smallest_eigenvalue    # most negative eigenvalue of sym(J_R)
result.largest_eigenvalue     # most positive eigenvalue of sym(J_R)
result.nearest_zero_eigenvalue  # eigenvalue closest to zero (stability margin)

# Print a diagnostic report:
println(PowerFlows.get_definiteness_report(jacobian, time_step))
```

### Monitoring During Solver Iterations

Enable per-iteration monitoring with the `monitor_jacobian` solver setting:

```julia
pf = ACPowerFlow{NewtonRaphsonACPowerFlow}(;
    solver_settings = Dict{Symbol, Any}(:monitor_jacobian => true),
)
result = solve_power_flow(pf, sys)
```
