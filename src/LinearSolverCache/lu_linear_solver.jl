mutable struct LULinSolveCache <: LinearSolverCache
    lu_fact::SparseArrays.UMFPACK.UmfpackLU{Float64, Int32}
    reuse_symbolic::Bool
    check_pattern::Bool
end

LULinSolveCache(A::SparseMatrixCSC{Float64, Int32},
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true) =
    LULinSolveCache(LinearAlgebra.lu(A), reuse_symbolic, check_pattern)

function numeric_refactor!(cache::LULinSolveCache, A::SparseMatrixCSC{Float64, Int32})
    LinearAlgebra.lu!(cache.lu_fact, A; check = cache.check_pattern,
        reuse_symbolic = cache.reuse_symbolic)
end

function symbolic_refactor!(cache::LULinSolveCache, A::SparseMatrixCSC{Float64, Int32})
    LinearAlgebra.lu!(cache.lu_fact, A; reuse_symbolic = false)
end

# if subtype doesn't define symbolic_factor!, default to calling symbolic_refactor!.
symbolic_factor!(cache::LULinSolveCache, A::SparseMatrixCSC{Float64, Int32}) =
    symbolic_refactor!(cache, A)

function solve!(cache::LULinSolveCache, B::Vector{Float64})
    LinearAlgebra.ldiv!(cache.lu_fact, B)
end
