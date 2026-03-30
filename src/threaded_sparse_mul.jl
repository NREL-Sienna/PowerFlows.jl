"""
    THREADED_MUL_MIN_DIM

Minimum matrix dimension (number of columns for transpose, rows for standard) below which
`threaded_mul!` falls back to `LinearAlgebra.mul!`. Threading overhead dominates for small
problems, so this threshold avoids slowdowns on systems with fewer buses/arcs.
"""
const THREADED_MUL_MIN_DIM = 1000

"""
    threaded_mul!(y, transpose(A), x) -> y
    threaded_mul!(y, adjoint(A), x)  -> y

Compute `y = Aᵀx` (or `y = A'x`) using multi-threaded parallelism for
`A::SparseMatrixCSC`.

In CSC (Compressed Sparse Column) format, computing the **transpose**-vector product is
naturally data-parallel: each output element `y[i]` depends only on column `i` of `A`,
so threads never contend on writes. This is in contrast to the standard product `y = Ax`,
where multiple columns scatter into overlapping rows of `y`.

Falls back to single-threaded `LinearAlgebra.mul!` when:
- the output dimension is below `THREADED_MUL_MIN_DIM`, or
- only one Julia thread is available (`Threads.nthreads() <= 1`).

!!! note
    The proposed function in the original issue computes `Aᵀx`, **not** `Ax`, despite
    the generic name `threaded_mul!`. It also unnecessarily requires a square matrix.
    This implementation lifts that restriction and dispatches on `Transpose`/`Adjoint`
    wrappers for type safety.
"""
function threaded_mul!(
    y::AbstractVector,
    At::LinearAlgebra.Transpose{<:Real, <:SparseMatrixCSC},
    x::AbstractVector,
)
    A = parent(At)
    n = A.n  # columns of A = length of y
    m = A.m  # rows of A = length of x

    if n < THREADED_MUL_MIN_DIM || Threads.nthreads() <= 1
        mul!(y, At, x)
        return y
    end

    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    Threads.@threads for i in 1:n
        tmp = zero(eltype(y))
        @inbounds for j in colptr[i]:(colptr[i + 1] - 1)
            tmp += nzval[j] * x[rowval[j]]
        end
        @inbounds y[i] = tmp
    end
    return y
end

function threaded_mul!(
    y::AbstractVector,
    At::LinearAlgebra.Adjoint{<:Real, <:SparseMatrixCSC},
    x::AbstractVector,
)
    # For real matrices, adjoint and transpose are identical.
    return threaded_mul!(y, transpose(parent(At)), x)
end

"""
    threaded_mul!(Y, transpose(A), X) -> Y

Multi-threaded transpose sparse-matrix × dense-matrix product `Y = AᵀX`.

Parallelises over rows of `Y` (i.e., columns of `A`). Each thread computes a full row
of `Y` independently, so there are no write conflicts.
"""
function threaded_mul!(
    Y::AbstractMatrix,
    At::LinearAlgebra.Transpose{<:Real, <:SparseMatrixCSC},
    X::AbstractMatrix,
)
    A = parent(At)
    n = A.n           # output rows
    n_cols = size(X, 2)  # time steps / columns

    if n < THREADED_MUL_MIN_DIM || Threads.nthreads() <= 1
        mul!(Y, At, X)
        return Y
    end

    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    Threads.@threads for i in 1:n
        @inbounds for t in 1:n_cols
            tmp = zero(eltype(Y))
            for j in colptr[i]:(colptr[i + 1] - 1)
                tmp += nzval[j] * X[rowval[j], t]
            end
            Y[i, t] = tmp
        end
    end
    return Y
end

function threaded_mul!(
    Y::AbstractMatrix,
    At::LinearAlgebra.Adjoint{<:Real, <:SparseMatrixCSC},
    X::AbstractMatrix,
)
    return threaded_mul!(Y, transpose(parent(At)), X)
end

# ── Generic fallbacks for dense or other matrix types ──────────────────────
# When the wrapped matrix is not a SparseMatrixCSC (e.g. dense PTDF), delegate
# to LinearAlgebra.mul! which already uses optimised BLAS routines.

function threaded_mul!(
    y::AbstractVecOrMat,
    At::Union{LinearAlgebra.Transpose, LinearAlgebra.Adjoint},
    x::AbstractVecOrMat,
)
    mul!(y, At, x)
    return y
end
