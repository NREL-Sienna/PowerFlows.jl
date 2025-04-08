"""A cached linear solver using KLU. 
# Fields:
- `K`: the underlying KLU object.
- `reuse_symbolic::Bool`: reuse the symbolic factorization. Defaults to true.
- `check_pattern::Bool`: if true, `numeric_refactor!` verifies that the new
matrix has the same sparsity structure. Defaults to true.
- `been_factored::Bool`: used internally. Keeps track of whether we've done a numeric factorization."""
mutable struct KLULinSolveCache{T} <: LinearSolverCache{T}
    K::KLU.KLUFactorization{Float64, T}
    reuse_symbolic::Bool
    check_pattern::Bool
    been_factored::Bool
end

function KLULinSolveCache(
    A::SparseMatrixCSC{Float64, T},
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true,
) where {T <: TIs}
    return KLULinSolveCache(KLU.KLUFactorization(A), reuse_symbolic, check_pattern, false)
end

get_reuse_symbolic(cache::KLULinSolveCache) = cache.reuse_symbolic

"""Compares sparsity structure of A with that of cache. Note the fields of `cache.K` are 
0-indexed, while those of A are 1-indexed."""
function same_pattern(cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    # KLU does this same increment-compare-decrement pattern: it avoids allocations,
    # and ensures that the matrix is in a valid state when we return.
    KLU.increment!(cache.K.rowval)
    KLU.increment!(cache.K.colptr)
    res = A.colptr == cache.K.colptr && A.rowval == cache.K.rowval
    KLU.decrement!(cache.K.rowval)
    KLU.decrement!(cache.K.colptr)
    return res
end

"""Frees up the current symbolic and numeric factorizations stored by `cache`, if non-null.
Then computes the symbolic factorization of `A` and stores that to `cache`."""
function symbolic_factor!(
    cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    if !(size(A, 1) == cache.K.n && size(A, 2) == cache.K.n)
        throw(
            DimensionMismatch(
                "Can't factor: matrix has different dimensions."),
        )
    end
    if cache.check_pattern && !same_pattern(cache, A)
        cache.K = KLU.KLUFactorization(A)
    end
    KLU.klu_analyze!(cache.K)
    cache.been_factored = false
    return
end

"""Symbolic refactor. Behavior depends on the values of `cache.reuse_symbol` and 
`cache.check_pattern`. There are 3 cases:
- `!reuse_symbol`: always refactor. Just calls `symbolic_factor(cache, A)`.
- `reuse_symbol && check_pattern`: checks if the symbolic structure of `A` matches the
    cached one, and throws an error if it doesn't. This is to prevent bad input: we expected 
    the structure to be the same, but it isn't.
- `reuse_symbol && !check pattern`: do nothing. Assume the structure of `A` matches the cached
     one."""
function symbolic_refactor!(
    cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    if cache.reuse_symbolic && cache.check_pattern
        if !(size(A, 1) == cache.K.n && size(A, 2) == cache.K.n)
            throw(
                DimensionMismatch(
                    "Can't refactor: new matrix has different dimensions."),
            )
        end
        if !same_pattern(cache, A)
            throw(
                ArgumentError(
                    "Matrix has different sparse structure. Either make cache with reuse_" *
                    "symbolic = false, or call symbolic_factor! instead of symbolic_refactor!",
                ),
            )
        end
    elseif !cache.reuse_symbolic
        symbolic_factor!(cache, A)
    end
    return
end

"""Frees numeric factorization stored by `cache`, if non-null. If `cache.check_pattern` 
is `true` and the sparse matrix structure of `A` doesn't match the cached one, 
throws an error. Finally, computes the numeric factorization of `A` and stores that to 
`cache`."""
function numeric_refactor!(
    cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    # KLU differentiates between factor and refactor. Might be a better solution
    # via checking the status of the KLU object.
    if !cache.been_factored
        KLU.klu_factor!(cache.K)
        cache.been_factored = true
    end
    # this checks that the symbolic structure matches for us.
    KLU.klu!(cache.K, A)
    return
end

"""Solves Ax = B in-place, where B is a vector or matrix."""
function solve!(
    cache::KLULinSolveCache{T},
    B::StridedVecOrMat{Float64},
) where {T <: TIs}
    # this checks the dimensions for us.
    return KLU.solve!(cache.K, B)
end

"""Solve with iterative refinement."""
function solve_w_refinement(cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
    B::StridedVecOrMat{Float64},
    tol::Float64 = 1e-6,
) where {T <: TIs}
    # parameters and strategy here could be fine-tuned.
    bNorm = norm(B, 1)
    XB = zeros(size(B))
    r = B - A * XB
    MAX_ITERS = DEFAULT_REFINEMENT_MAX_ITER
    iters = 0
    while iters < MAX_ITERS && norm(r, 1) >= bNorm * tol #*cache.condest
        lastError = norm(r, 1)
        solve!(cache, r)
        XB .+= r
        r .= B - A * XB
        iters += 1
        if norm(r, 1) > lastError
            @error("Iterative refinement failed: error is getting worse.")
            return XB
        end
    end
    @debug("Iterative refined converged in $iters iterations.")
    return XB
end
