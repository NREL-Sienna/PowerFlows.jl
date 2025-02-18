mutable struct LULinSolveCache <: LinearSolverCache
    lu_fact::SparseArrays.UMFPACK.UmfpackLU{Float64, Int32}
    reuse_symbolic::Bool
    check_pattern::Bool
end

function LULinSolveCache(A::SparseMatrixCSC{Float64, Int32},
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true)
    LULinSolveCache(LinearAlgebra.lu(A), reuse_symbolic, check_pattern)
end

function numeric_refactor!(
    cache::LULinSolveCache,
    A::SparseMatrixCSC{Float64, Int32},
)
    LinearAlgebra.lu!(cache.lu_fact, A; check = cache.check_pattern,
        reuse_symbolic = cache.reuse_symbolic)
end

function symbolic_refactor!(
    cache::LULinSolveCache,
    A::SparseMatrixCSC{Float64, Int32},
) 
    LinearAlgebra.lu!(cache.lu_fact, A; reuse_symbolic = false)
end

function symbolic_factor!(
    cache::LULinSolveCache,
    A::SparseMatrixCSC{Float64, Int32},
) 
    symbolic_refactor!(cache, A)
end

function solve!(cache::LULinSolveCache, B::StridedVecOrMat{Float64})
    LinearAlgebra.ldiv!(cache.lu_fact, B)
end
