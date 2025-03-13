# select correct function from KLU based on type.
__free_sym(::Type{Int32}, args...) = KLU.klu_free_symbolic(args...)
__free_sym(::Type{Int64}, args...) = KLU.klu_l_free_symbolic(args...)

__free_num(::Type{Int32}, args...) = KLU.klu_free_numeric(args...)
__free_num(::Type{Int64}, args...) = KLU.klu_l_free_numeric(args...)

__analyze(::Type{Int32}, args...) = KLU.klu_analyze(args...)
__analyze(::Type{Int64}, args...) = KLU.klu_l_analyze(args...)

__factor(::Type{Int32}, args...) = KLU.klu_factor(args...)
__factor(::Type{Int64}, args...) = KLU.klu_l_factor(args...)

__refactor(::Type{Int32}, args...) = KLU.klu_refactor(args...)
__refactor(::Type{Int64}, args...) = KLU.klu_l_refactor(args...)

__solve(::Type{Int32}, args...) = KLU.klu_solve(args...)
__solve(::Type{Int64}, args...) = KLU.klu_l_solve(args...)

__symType(::Type{Int32}) = KLU.klu_symbolic
__symType(::Type{Int64}) = KLU.klu_l_symbolic
__numType(::Type{Int32}) = KLU.klu_numeric
__numType(::Type{Int64}) = KLU.klu_l_numeric

"""A cached linear solver using KLU. Carefully written so as to minimize
allocations: solve! and numeric_refactor! are completely non-allocating.
# Fields:
- `K`: the underlying KLU object.
- `reuse_symbolic::Bool`: reuse the symbolic factorization. Defaults to true.
- `check_pattern::Bool`: if true, `numeric_refactor!` verifies that the new
matrix has the same sparsity structure. Defaults to true.
-`rf_common`, `rf_symbolic`, `rf_numeric`: internal usage. Stored to avoid allocations."""
mutable struct KLULinSolveCache{T} <: LinearSolverCache{T}
    K::KLU.KLUFactorization{Float64, T}
    reuse_symbolic::Bool
    check_pattern::Bool
    # condest::Float64 # condition number. could be useful for interative refinement.
    rf_common::Union{Base.RefValue{KLU.klu_common}, Base.RefValue{KLU.klu_l_common}}
    # klu_free_{numeric/symbolic} requires these. store them to avoid allocating.
    rf_symbolic::Union{
        Base.RefValue{Ptr{KLU.klu_symbolic}},
        Base.RefValue{Ptr{KLU.klu_l_symbolic}},
    }
    rf_numeric::Union{
        Base.RefValue{Ptr{KLU.klu_numeric}},
        Base.RefValue{Ptr{KLU.klu_l_numeric}},
    }
end

"""Constructor."""
function KLULinSolveCache(
    A::SparseMatrixCSC{Float64, T},
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true,
) where {T <: TIs}
    symType = __symType(T)
    numType = __numType(T)
    K = KLU.KLUFactorization(A)
    return KLULinSolveCache(K, reuse_symbolic, check_pattern,
        # NaN,
        Ref(K.common), Ref(Ptr{symType}(K._symbolic)),
        Ref(Ptr{numType}(K._numeric)))
end

get_reuse_symbolic(cache::KLULinSolveCache) = cache.reuse_symbolic

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
    # cache.condest = NaN # will need to re-calculate condition number.
    if cache.K._symbolic != C_NULL
        @assert cache.rf_symbolic != Ref(C_NULL)
        __free_sym(T, cache.rf_symbolic, cache.rf_common)
        cache.K._symbolic = C_NULL
    end
    if cache.K._numeric != C_NULL
        @assert cache.rf_numeric != Ref(C_NULL)
        __free_num(T, cache.rf_numeric, cache.rf_common)
        cache.K._numeric = C_NULL
    end

    # copy over new info in a minimally-allocating way.
    resize!(cache.K.colptr, length(A.colptr))
    copyto!(cache.K.colptr, A.colptr)
    KLU.decrement!(cache.K.colptr)

    resize!(cache.K.rowval, length(A.rowval))
    copyto!(cache.K.rowval, A.rowval)
    KLU.decrement!(cache.K.rowval)

    resize!(cache.K.nzval, length(A.nzval))
    @assert cache.K._symbolic == C_NULL

    sym = __analyze(T, cache.K.n, cache.K.colptr, cache.K.rowval, cache.rf_common)
    if sym != C_NULL
        symType = __symType(T)
        cache.K._symbolic = sym
        cache.rf_symbolic[] = Ptr{symType}(cache.K._symbolic)
    end
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
        # KLU does this same increment-compare-decrement pattern.
        KLU.increment!(cache.K.rowval)
        KLU.increment!(cache.K.colptr)
        shouldErr::Bool = A.colptr != cache.K.colptr || A.rowval != cache.K.rowval
        KLU.decrement!(cache.K.rowval)
        KLU.decrement!(cache.K.colptr)
        shouldErr && throw(
            ArgumentError(
                "Matrix has different sparse structure. Either make cache with reuse_" *
                "symbolic = false, or call symbolic_factor! instead of symbolic_refactor!",
            ),
        )
    elseif !cache.reuse_symbolic
        symbolic_factor!(cache, A) # this also sets rf_symbolic.
        @assert cache.rf_symbolic[] == cache.K._symbolic
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
    if cache.K._symbolic == C_NULL
        @error("Need to call symbolic_factor! first.")
    end
    # cache.condest = NaN # will need to re-calculate condition number.
    if cache.K._numeric == C_NULL
        cache.K._numeric = __factor(T,
            cache.K.colptr,
            cache.K.rowval,
            A.nzval,
            cache.K._symbolic,
            cache.rf_common)
        if cache.K._numeric == C_NULL
            @warn("factor failed")
            KLU.kluerror(cache.K.common)
        end
    else
        if cache.check_pattern
            # 0 vs 1-indexing. If mismatch, put things back before throwing error.
            KLU.increment!(cache.K.rowval)
            KLU.increment!(cache.K.colptr)
            shouldErr::Bool = A.colptr != cache.K.colptr || A.rowval != cache.K.rowval
            KLU.decrement!(cache.K.rowval)
            KLU.decrement!(cache.K.colptr)
            shouldErr && throw(
                ArgumentError(
                    "Cannot numeric_refactor: matrix has different sparse structure.",
                ),
            )
        end
        ok = __refactor(T,
            cache.K.colptr,
            cache.K.rowval,
            A.nzval,
            cache.K._symbolic,
            cache.K._numeric,
            cache.rf_common,
        )
        if (ok != 1)
            @warn("refactor failed")
            KLU.kluerror(cache.K.common)
        end
        numType = __numType(T)
        cache.rf_numeric[] = Ptr{numType}(cache.K._numeric)
    end
    return
end

function solve!(cache::KLULinSolveCache{T}, B::StridedVecOrMat{Float64}) where {T <: TIs}
    size(B, 1) == cache.K.n || throw(
        DimensionMismatch(
            "Need size(B, 1) to equal $(cache.K.n), but got $(size(B, 1))."),
    )
    stride(B, 1) == 1 || throw(ArgumentError("B must have unit strides"))
    isok = __solve(T, cache.K._symbolic, cache.K._numeric,
        size(B, 1), size(B, 2), B, cache.rf_common)
    isok == 0 && KLU.kluerror(cache.K.common)
    return B
end

function solve_w_refinement(cache::KLULinSolveCache{T}, A::SparseMatrixCSC{Float64, T},
    B::StridedVecOrMat{Float64}, tol::Float64 = 1e-6) where {T <: TIs}
    cache.K._numeric == C_NULL && @error("not factored yet")
    # update condition number, if needed
    #=if isnan(cache.condest)
        condestFcn = T <: Int32 ? KLU.klu_condest : KLU.klu_l_condest
        ok = condestFcn(cache.K.colptr, cache.K.nzval, cache.K._symbolic, cache.K._numeric, cache.rf_common)
        if ok == 0
            KLU.kluerror(cache.K.common)
        end
        cache.condest = K.common.condest
    end=#
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
