struct CholeskyHessianSolver <: HessianSolver
    F::SparseArrays.CHOLMOD.Factor{Float64, Int32}
    mat::FixedStructureCHOLMOD
    buff::Vector{Float64} # buffer for solving
end

function CholeskyHessianSolver(H::SparseMatrixCSC{Float64, Int32})
    mat = FixedStructureCHOLMOD(H)
    # I need to create a CHOLMOD factorization object, so I also symbolic factor it here.
    n = size(H, 1)
    return CholeskyHessianSolver(symbolic_factor(mat), mat, zeros(n))
end

function symbolic_factor!(::CholeskyHessianSolver, ::SparseMatrixCSC{Float64, Int32})
    # hSolver.F = symbolic_factor(hSolver.mat) # if I make FixedStructureCHOLMOD mutable.
    return
end

function modify_and_numeric_factor!(
    hSolver::CholeskyHessianSolver,
    H::SparseMatrixCSC{Float64, Int32},
)
    minDiagElem = minimum(H[i, i] for i in axes(H, 1))
    τ_old = 0.0
    if minDiagElem > 0.0
        τ = 0.0
    else
        τ = -minDiagElem + β
    end
    @debug "initial τ = $τ"
    nonsingular = false
    fill!(hSolver.buff, 1.0)
    while !nonsingular
        for i in axes(H, 1)
            H[i, i] += τ - τ_old # now try H + τ*I
        end
        try
            fill!(hSolver.buff, 1.0) # reset b to a vector of ones.
            set_values!(hSolver.mat, SparseArrays.nonzeros(H))
            numeric_factor!(hSolver.F, hSolver.mat)
            hSolver.F \ hSolver.buff # sometimes the error isn't thrown until we use the factorization.
            nonsingular = true
            @debug "nonsingular with τ = $τ"
        catch e
            if e isa SparseArrays.CHOLMOD.PosDefException
                τ_old = τ
                τ *= 2.0
                τ = max(τ, β) # ensure τ is at least β
            else
                rethrow(e)
            end
        end
    end
    return
end

function solve!(solver::CholeskyHessianSolver, b::Vector{Float64})
    copyto!(solver.buff, b)
    b .= solver.F \ solver.buff
end
