abstract type LinearSolverCache end

function symbolic_factor!(cache::LinearSolverCache, A::SparseMatrixCSC{Float64, Int64})
    throw(AbstractMethodError(:symbolic_factor!))
end

function symbolic_refactor!(cache::LinearSolverCache, A::SparseMatrixCSC{Float64, Int64})
    throw(AbstractMethodError(:symbolic_refactor!))
end

function numeric_refactor!(cache::LinearSolverCache,  A::SparseMatrixCSC{Float64, Int64})
    throw(AbstractMethodError(:numeric_refactor!))
end

function solve!(cache::LinearSolverCache, B::Vector{Float64})
    throw(AbstractMethodError(:solve!))
end

function full_refactor!(cache::LinearSolverCache, A::SparseMatrixCSC{Float64, Int64})
    symbolic_refactor!(cache, A)
    numeric_refactor!(cache, A)
    return
end

function full_factor!(cache::LinearSolverCache, A::SparseMatrixCSC{Float64, Int64})
    symbolic_factor!(cache, A)
    numeric_refactor!(cache, A)
    return
end