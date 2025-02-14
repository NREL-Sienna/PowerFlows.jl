const TIs = Union{Int32, Int64}
abstract type LinearSolverCache{T <: TIs} end
function symbolic_factor!(
    cache::LinearSolverCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    throw(AbstractMethodError(:symbolic_factor!))
end

function symbolic_refactor!(
    cache::LinearSolverCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    throw(AbstractMethodError(:symbolic_refactor!))
end

function numeric_refactor!(
    cache::LinearSolverCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    throw(AbstractMethodError(:numeric_refactor!))
end

function solve!(cache::LinearSolverCache{T}, B::Vector{Float64}) where {T <: TIs}
    throw(AbstractMethodError(:solve!))
end

function full_refactor!(
    cache::LinearSolverCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    symbolic_refactor!(cache, A)
    numeric_refactor!(cache, A)
    return
end

function full_factor!(
    cache::LinearSolverCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: TIs}
    symbolic_factor!(cache, A)
    numeric_refactor!(cache, A)
    return
end
