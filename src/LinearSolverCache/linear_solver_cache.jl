const TIs = Union{Int32, Int64}
"""Abstract supertype for all cached linear solvers.
Subtypes must implement: `symbolic_factor!`, `symbolic_refactor!`,
`numeric_refactor!` (which doubles as `numeric_factor!`), and `solve!`."""
abstract type LinearSolverCache{T <: TIs} end

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
