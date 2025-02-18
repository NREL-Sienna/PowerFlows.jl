mutable struct LULinSolveCache{T} <: LinearSolverCache{T}
    lu_fact::SparseArrays.UMFPACK.UmfpackLU{Float64, T}
    reuse_symbolic::Bool
    check_pattern::Bool
end

function LULinSolveCache(A::SparseMatrixCSC{Float64, T},
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true) where {T <: TIs}
    LULinSolveCache(LinearAlgebra.lu(A), reuse_symbolic, check_pattern)
end

function numeric_refactor!(
    cache::LULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    LinearAlgebra.lu!(cache.lu_fact, A; check = cache.check_pattern,
        reuse_symbolic = cache.reuse_symbolic)
end

function symbolic_refactor!(
    cache::LULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    LinearAlgebra.lu!(cache.lu_fact, A; reuse_symbolic = false)
end

function symbolic_factor!(
    cache::LULinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    symbolic_refactor!(cache, A)
end

function solve!(cache::LULinSolveCache{T}, B::StridedVecOrMat{Float64}) where {T <: TIs}
    LinearAlgebra.ldiv!(cache.lu_fact, B)
end
