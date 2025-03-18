"""Abstract supertype for all cached linear solvers.
Subtypes must implement: `symbolic_factor!`, `symbolic_refactor!`,
`numeric_refactor!` (which doubles as `numeric_factor!`), and `solve!`."""
abstract type LinearSolverCache{T <: Integer} end

function full_refactor!(
    cache::LinearSolverCache,
    A::SparseMatrixCSC{Float64},
)
    symbolic_refactor!(cache, A)
    numeric_refactor!(cache, A)
    return
end

function full_factor!(
    cache::LinearSolverCache,
    A::SparseMatrixCSC{Float64},
)
    symbolic_factor!(cache, A)
    numeric_refactor!(cache, A)
    return
end
