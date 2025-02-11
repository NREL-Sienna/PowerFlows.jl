using LinearAlgebra

mutable struct KLULinSolveCache <: LinearSolverCache
    K::KLU.KLUFactorization{Float64, Int64}
    reuse_symbolic::Bool
    check_pattern::Bool
    rf_common::Base.RefValue{KLU.klu_l_common}
    # klu_l_free_{numeric/symbolic} requires these. store them to avoid allocating.
    rf_symbolic::Base.RefValue{Ptr{KLU.klu_l_symbolic}}
    rf_numeric::Base.RefValue{Ptr{KLU.klu_l_numeric}}
end

function KLULinSolveCache(
    A::SparseMatrixCSC{Float64, Int64},
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true
)
    K = KLU.KLUFactorization(A)
    # KLU.kluerror(K.common) # necessary?
    return KLULinSolveCache(K, reuse_symbolic, check_pattern,
                        Ref(K.common), Ref(Ptr{KLU.klu_l_symbolic}(K._symbolic)),
                                        Ref(Ptr{KLU.klu_l_numeric}(K._numeric)))
end

get_reuse_symbolic(cache::KLULinSolveCache) = cache.reuse_symbolic

# only does the symbolic, not the numeric.
function symbolic_factor!(cache::KLULinSolveCache, A::SparseMatrixCSC{Float64, Int64})
    @assert size(A, 1) == cache.K.n && size(A, 2) == cache.K.n

    if cache.K._symbolic != C_NULL
        @assert cache.rf_symbolic != Ref(C_NULL)
        KLU.klu_l_free_symbolic(cache.rf_symbolic, cache.rf_common)
        cache.K._symbolic = C_NULL
        @assert cache.rf_symbolic[] == C_NULL
    end
    if cache.K._numeric != C_NULL
        @assert cache.rf_numeric != Ref(C_NULL)
        KLU.klu_l_free_numeric(cache.rf_numeric, cache.rf_common)
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
    sym = KLU.klu_l_analyze(cache.K.n, cache.K.colptr, cache.K.rowval, cache.rf_common)
    if sym != C_NULL
        cache.K._symbolic = sym
        cache.rf_symbolic[] = Ptr{KLU.klu_l_symbolic}(cache.K._symbolic)
        # we do need the above line: without it, this assert triggers
        @assert cache.rf_symbolic[] == cache.K._symbolic
    end
    return
end

function symbolic_refactor!(cache::KLULinSolveCache, A::SparseMatrixCSC{Float64, Int64})
    # case 1: reuse_symbol && !check_pattern => assume pattern hasn't changed.
    # case 2: reuse_symbol && check_pattern => check pattern
    # case 3: !reuse_symbol => refactor.
    if cache.reuse_symbolic && cache.check_pattern
        @assert size(A, 1) == cache.K.n && size(A, 2) == cache.K.n
        # KLU does this same increment-compare-decrement pattern.
        KLU.increment!(cache.K.rowval)
        KLU.increment!(cache.K.colptr)
        shouldErr::Bool =  A.colptr != cache.K.colptr || A.rowval != cache.K.rowval
        KLU.decrement!(cache.K.rowval)
        KLU.decrement!(cache.K.colptr)
        shouldErr && throw(ArgumentError("Matrix has different sparse structure. " *
                                "Either make cache with reuse_symbolic = false, " *
                        "or call symbolic_factor! [instead of symbolic_refactor!]."))
    elseif !cache.reuse_symbolic
        symbolic_factor!(cache, A) # this also sets rf_symbolic.
        @assert cache.rf_symbolic[] == cache.K._symbolic
    end
    return
end

function numeric_refactor!(cache::KLULinSolveCache, A::SparseMatrixCSC{Float64, Int64})

    if cache.K._symbolic == C_NULL
        @error("Need to call symbolic_refactor! first.")
    end
    if cache.K._numeric == C_NULL
        cache.K._numeric = KLU.klu_l_factor(cache.K.colptr,
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
            shouldErr::Bool =  A.colptr != cache.K.colptr || A.rowval != cache.K.rowval
            KLU.decrement!(cache.K.rowval)
            KLU.decrement!(cache.K.colptr)
            shouldErr && throw(ArgumentError("Cannot numeric_refactor: "*
                                        "matrix has different sparse structure."))
        end
        ok = KLU.klu_l_refactor(cache.K.colptr,
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
        cache.rf_numeric[] = Ptr{KLU.klu_l_numeric}(cache.K._numeric)
        # we do need the above line: without it, the below assert triggers.
        @assert cache.rf_numeric[] == Ptr{KLU.klu_l_numeric}(cache.K._numeric)
    end
    return
end

function solve!(cache::KLULinSolveCache, B::Vector{Float64})
    if length(B) != cache.K.n
        throw(LinearAlgebra.DimensionMismatch("b must be a"*
                                    " vector of length $(cache.K.n)."))
    end
    isok = KLU.klu_l_solve(cache.K._symbolic, cache.K._numeric,
                            size(B, 1), size(B, 2), B, cache.rf_common)
    isok == 0 && KLU.kluerror(cache.K.common)
    return B
end
