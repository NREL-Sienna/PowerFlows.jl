using LinearAlgebra
mutable struct KLULinSolveCache{T} <: LinearSolverCache{T}
    K::KLU.KLUFactorization{Float64, T}
    reuse_symbolic::Bool
    check_pattern::Bool
    # condest::Float64
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

function KLULinSolveCache(
    A::SparseMatrixCSC{Float64, T},
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true,
) where {T <: TIs}
    symType = T == Int32 ? KLU.klu_symbolic : KLU.klu_l_symbolic
    numType = T == Int32 ? KLU.klu_numeric : KLU.klu_l_numeric
    K = KLU.KLUFactorization(A)
    # KLU.kluerror(K.common) # necessary?
    return KLULinSolveCache(K, reuse_symbolic, check_pattern,
        # NaN,
        Ref(K.common), Ref(Ptr{symType}(K._symbolic)),
        Ref(Ptr{numType}(K._numeric)))
end

get_reuse_symbolic(cache::KLULinSolveCache{T}) where {T <: TIs} = cache.reuse_symbolic

# only does the symbolic, not the numeric.
function symbolic_factor!(
    cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    # cache.condest = NaN # will need to re-calculate condition number.
    @assert size(A, 1) == cache.K.n && size(A, 2) == cache.K.n

    if cache.K._symbolic != C_NULL
        freeFcn = T == Int32 ? KLU.klu_free_symbolic : KLU.klu_l_free_symbolic
        @assert cache.rf_symbolic != Ref(C_NULL)
        freeFcn(cache.rf_symbolic, cache.rf_common)
        cache.K._symbolic = C_NULL
        @assert cache.rf_symbolic[] == C_NULL
    end
    if cache.K._numeric != C_NULL
        freeFcn = T == Int32 ? KLU.klu_free_numeric : KLU.klu_l_free_numeric
        @assert cache.rf_numeric != Ref(C_NULL)
        freeFcn(cache.rf_numeric, cache.rf_common)
        cache.K._numeric = C_NULL
        @assert cache.rf_numeric[] == C_NULL
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

    analyzeFcn = T == Int32 ? KLU.klu_analyze : KLU.klu_l_analyze
    sym = analyzeFcn(cache.K.n, cache.K.colptr, cache.K.rowval, cache.rf_common)
    if sym != C_NULL
        symType = T == Int32 ? KLU.klu_symbolic : KLU.klu_l_symbolic
        cache.K._symbolic = sym
        cache.rf_symbolic[] = Ptr{symType}(cache.K._symbolic)
        # we do need the above line: without it, this assert triggers
        @assert cache.rf_symbolic[] == cache.K._symbolic
    end
    return
end

function symbolic_refactor!(
    cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    # case 1: reuse_symbol && !check_pattern => assume pattern hasn't changed.
    # case 2: reuse_symbol && check_pattern => check pattern
    # case 3: !reuse_symbol => refactor.
    if cache.reuse_symbolic && cache.check_pattern
        @assert size(A, 1) == cache.K.n && size(A, 2) == cache.K.n
        # KLU does this same increment-compare-decrement pattern.
        KLU.increment!(cache.K.rowval)
        KLU.increment!(cache.K.colptr)
        shouldErr::Bool = A.colptr != cache.K.colptr || A.rowval != cache.K.rowval
        KLU.decrement!(cache.K.rowval)
        KLU.decrement!(cache.K.colptr)
        shouldErr && throw(
            ArgumentError(
                "Matrix has different sparse structure. " *
                "Either make cache with reuse_symbolic = false, " *
                "or call symbolic_factor! [instead of symbolic_refactor!].",
            ),
        )
    elseif !cache.reuse_symbolic
        symbolic_factor!(cache, A) # this also sets rf_symbolic.
        @assert cache.rf_symbolic[] == cache.K._symbolic
    end
    return
end

function numeric_refactor!(
    cache::KLULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    if cache.K._symbolic == C_NULL
        @error("Need to call symbolic_factor! first.")
    end
    # cache.condest = NaN # will need to re-calculate condition number.
    if cache.K._numeric == C_NULL
        factorFcn = T == Int32 ? KLU.klu_factor : KLU.klu_l_factor
        cache.K._numeric = factorFcn(cache.K.colptr,
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
                    "Cannot numeric_refactor: " *
                    "matrix has different sparse structure.",
                ),
            )
        end
        refactorFcn = T == Int32 ? KLU.klu_refactor : KLU.klu_l_refactor
        ok = refactorFcn(cache.K.colptr,
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
        numType = T == Int32 ? KLU.klu_numeric : KLU.klu_l_numeric
        cache.rf_numeric[] = Ptr{numType}(cache.K._numeric)
        # we do need the above line: without it, the below assert triggers.
        @assert cache.rf_numeric[] == Ptr{numType}(cache.K._numeric)
    end
    return
end

function solve!(cache::KLULinSolveCache{T}, B::StridedVecOrMat{Float64}) where {T <: TIs}
    size(B, 1) == cache.K.n || throw(
        LinearAlgebra.DimensionMismatch(
            "Need size(B, 1) to equal $(cache.K.n), but got $(size(B, 1))."),
    )
    stride(B, 1) == 1 || throw(ArgumentError("B must have unit strides"))
    solveFcn = T == Int32 ? KLU.klu_solve : KLU.klu_l_solve
    isok = solveFcn(cache.K._symbolic, cache.K._numeric,
        size(B, 1), size(B, 2), B, cache.rf_common)
    isok == 0 && KLU.kluerror(cache.K.common)
    return B
end

function solve_w_refinement(cache::KLULinSolveCache{T}, A::SparseMatrixCSC{Float64, T},
    B::StridedVecOrMat{Float64}, tol::Float64 = 1e-6) where {T <: TIs}
    cache.K._numeric == C_NULL && @error("not factored yet")
    # update condition number, if needed
    #=if isnan(cache.condest)
        condestFcn = T == Int32 ? KLU.klu_condest : KLU.klu_l_condest
        ok = condestFcn(cache.K.colptr, cache.K.nzval, cache.K._symbolic, cache.K._numeric, cache.rf_common)
        if ok == 0
            KLU.kluerror(cache.K.common)
        end
        cache.condest = K.common.condest
    end=#
    bNorm = norm(B, 1)
    XB = zeros(size(B))
    r = B - A * XB
    MAX_ITERS = 10
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
