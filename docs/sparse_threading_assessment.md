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

## What This PR Implements

The proposed function's core idea (threaded transpose SpMV over CSC columns) is
sound but needed corrections. This PR takes the corrected concept and applies it
across the codebase:

### 1. `src/threaded_sparse_mul.jl` — Corrected `threaded_mul!`

- Dispatches on `Transpose`/`Adjoint` wrappers for type safety (not raw `SparseMatrixCSC`)
- No square-matrix restriction
- Size threshold (`THREADED_MUL_MIN_DIM = 1000`) to avoid overhead on small systems
- Automatic fallback to `LinearAlgebra.mul!` when single-threaded or below threshold
- Both vector (`Aᵀx`) and matrix (`AᵀX`) variants for multi-period support
- Generic fallback for dense matrices (delegates to BLAS)

### 2. DC power flow integration (`solve_dc_power_flow.jl`)

- `PTDFPowerFlowData`: transpose multiply now uses `threaded_mul!`
- `ABAPowerFlowData`: transpose multiply now uses `threaded_mul!`
- Loss factor computation: eliminates `Diagonal` allocation (addresses existing
  `PERF` comment), decomposes into fused 3-step computation

### 3. VirtualPTDF threading (`common.jl`)

- `my_mul_mt` vector and matrix variants now use `@threads` for large systems
  (gated on `THREADED_MUL_MIN_DIM`)

### 4. Multi-threaded Jacobian construction (`ac_power_flow_jacobian.jl`)

- `_update_jacobian_matrix_values!` uses `Base.Threads.@threads` over the
  per-bus loop, with a thread-local `MVector{4, Float64}` diagonal accumulator
  per bus (zero-cost stack allocation)
- Thread safety guaranteed: each bus writes exclusively to its own Jacobian
  rows (`2i-1`, `2i`)
- Gated on `THREADED_MUL_MIN_DIM` to avoid threading overhead on small systems
