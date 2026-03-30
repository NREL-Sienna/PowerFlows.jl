# Assessment: Threaded Sparse Matrix-Vector Multiplication for PowerFlows.jl

## Proposed Function

```julia
function threaded_mul!(y::Vector{T}, A::SparseMatrixCSC{T}, x::Vector{T}) where T <: Number
  A.m == A.n || error("A is not a square matrix!")
  @threads for i = 1 : A.n
    tmp = zero(T)
    @inbounds for j = A.colptr[i] : (A.colptr[i+1] - 1)
      tmp += A.nzval[j] * x[A.rowval[j]]
    end
    @inbounds y[i] = tmp
  end
  return y
end
```

## Critical Issues with the Proposed Code

### 1. Computes `Aᵀx`, not `Ax` (Bug)

The function iterates over **columns** of the CSC matrix and accumulates
`nzval[j] * x[rowval[j]]` into `y[i]`. In CSC format, column `i` stores
entries `A[rowval[j], i]`, so this computes:

```
y[i] = Σⱼ A[rowval[j], i] * x[rowval[j]] = (Aᵀx)[i]
```

This is **transpose multiplication**, not standard `Ax`. The function name
`threaded_mul!` is misleading. For non-symmetric matrices (e.g., the power flow
Jacobian), this gives incorrect results.

### 2. Unnecessary Square Matrix Restriction

The `A.m == A.n` check excludes rectangular matrices like the PTDF matrix
(`n_arcs × n_buses`), which is one of the main matrix-vector products in the
DC power flow path.

### 3. Thread Overhead vs. Problem Size

Power flow Jacobians are typically `2n × 2n` where `n` is the number of buses.
For systems under ~5,000 buses (10,000 rows), the `@threads` scheduling
overhead (~1–5 μs per spawn) can exceed the time saved. Julia's `@threads`
uses a static partitioning scheme that doesn't account for load imbalance from
varying column densities.

### 4. Missing Type Flexibility

The signature constrains `A`, `x`, and `y` to share the same element type `T`,
but `SparseMatrixCSC{Float64, Int32}` (this codebase uses `J_INDEX_TYPE = Int32`)
requires the index type parameter too. More importantly, several call sites use
`transpose(J.Jv)` wrappers, which are `Transpose{...}` types—not raw
`SparseMatrixCSC`.

## Where Sparse Matrix Multiplications Occur

| Call Site | Operation | Size | Frequency |
|-----------|-----------|------|-----------|
| `gradient_descent_ac_power_flow.jl:75` | `Jᵀ · F` (gradient) | `2n × 2n` | Every Adam iteration |
| `power_flow_method.jl:133` | `Jᵀ · r` (Cauchy point) | `2n × 2n` | Every trust region step |
| `power_flow_method.jl:267` | `J · Δx` (predicted residual) | `2n × 2n` | Every trust region step |
| `solve_dc_power_flow.jl:44` | `PTDFᵀ · P` (DC flows) | `arcs × buses` | Once per DC solve |
| `solve_dc_power_flow.jl:130` | `BAᵀ · θ` (ABA flows) | `arcs × buses` | Once per DC solve |
| `solve_dc_power_flow.jl:291` | `PTDF · Diag(R) · PTDFᵀ · P` (loss factors) | `arcs × buses` | Once if loss factors enabled |
| `common.jl:384` | `Xᵀ · row` (vPTDF multiply) | `buses × ts` | Per arc in vPTDF |

### Analysis

- **AC power flow (Newton/Trust Region):** The dominant cost is the KLU sparse
  factorization (`O(nnz^{1.5})` typically), not the SpMV (`O(nnz)`). Threading
  the SpMV would save a small fraction of total iteration time.
- **DC power flow:** The `PTDF' * P` multiplication is called once per solve,
  not in a loop. Threading has negligible impact.
- **Gradient descent:** The `Jᵀ · F` multiply is called every iteration and is
  the main linear algebra cost (no factorization). This is the best candidate
  for threading, but the Adam solver is a fallback method, not the primary path.

## Recommendation: Do Not Add the Proposed Function As-Is

The function has a correctness bug (computes transpose), an unnecessary
restriction (square only), and targets operations that are not the performance
bottleneck. The threading overhead will likely **hurt** performance for typical
power system sizes (< 10k buses).

### If Threading SpMV is Desired in the Future

For large-scale systems (50k+ buses), a correct implementation would:

1. **Use `polyester` or `@spawn`-based chunking** instead of `@threads` for
   lower overhead.
2. **Compute `Ax` correctly** by parallelizing over rows (requires CSR or
   a row-partitioned approach), or use the transpose formulation explicitly.
3. **Gate on problem size** — only thread when `nnz(A) > threshold`.
4. **Benchmark with realistic systems** — IEEE 118, Polish 2383, PEGASE 9241,
   and larger, to find the crossover point.

### What Would Actually Help Performance

Based on the codebase analysis, higher-impact optimizations would be:

1. **Multi-threaded Jacobian construction** (`ac_power_flow_jacobian.jl`):
   The per-bus neighbor loop that fills the Jacobian is embarrassingly parallel
   and dominates the non-factorization cost. **Implemented in this PR** — the
   `_update_jacobian_matrix_values!` function now uses `Base.Threads.@threads`
   over the bus loop, with a thread-local `MVector{4, Float64}` diagonal
   accumulator per bus (zero-cost stack allocation). Thread safety is guaranteed
   because each bus writes exclusively to its own Jacobian rows (`2i-1`, `2i`).
2. **Parallel multi-period solves** (`solve_dc_power_flow.jl`): Each time step
   is independent after factorization; threading across time steps would give
   near-linear speedup.
3. **Threaded vPTDF row computation** (`common.jl:362–387`): Each row of the
   VirtualPTDF is computed independently; `@threads` over arcs would help.
