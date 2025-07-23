@enum MUMPS_JOB begin
    INIT = -1
    CLEANUP = -2
    ANALYZE = 1
    FACTOR = 2
    SOLVE = 3
    SAVE = 7
    RESTORE = 8
    DELETE = -3
    FACTOR_CLEANUP = -4
end

function mumps_job!(mumps::Mumps, job::MUMPS_JOB)
    MUMPS.set_job!(mumps, Integer(job))
    MUMPS.invoke_mumps!(mumps)
    return
end

struct MUMPSHessianSolver <: HessianSolver
    mumps::Mumps
    function MUMPSHessianSolver(::SparseMatrixCSC)
        if !MPI.Initialized()
            MPI.Init()
        end
        icntl = deepcopy(MUMPS.default_icntl)
        icntl[4] = 1 # report errors only
        icntl[13] = 1 # due to our _modify_hessian! strategy,
        # we need to know the exact negative pivots.
        mumps = Mumps{Float64}(MUMPS.mumps_symmetric, icntl, MUMPS.default_cntl32)
        MPI.add_finalize_hook!(() -> MUMPS.finalize(mumps))
        return new(mumps)
    end
end

function symbolic_factor!(mumpsSolver::MUMPSHessianSolver, H::SparseMatrixCSC)
    MUMPS.associate_matrix!(mumpsSolver.mumps, H)
    mumps_job!(mumpsSolver.mumps, ANALYZE)
    if mumpsSolver.mumps.infog[1] < 0
        throw(ErrorException("MUMPS symbolic factorization failed with "*
                "error code $(mumpsSolver.mumps.infog[1])"))
    end
end


function modify_and_numeric_factor!(mumpsSolver::MUMPSHessianSolver, H::SparseMatrixCSC)
    minDiagElem = minimum(H[i, i] for i in axes(H, 1))
    if minDiagElem > 0.0
        τ = 0.0
    else
        τ = -minDiagElem + β
        for i in axes(H, 1)
            H[i, i] += τ
        end
    end
    @debug "initial τ = $τ"
    # cleanup, if necessary, before refactoring.
    if Int(mumpsSolver.mumps.job) in Int.([FACTOR, SOLVE])
        mumps_job!(mumpsSolver.mumps, FACTOR_CLEANUP)
    end
    # PERF: pass pointers to mumps, so we don't need to associate_matrix! each time.
    #       MUMPS stores things in COO format, though: just pass a pointer to 
    #       nzval? That's passing a pointer to a struct internals, though... 
    #       See issue #160 in the MUMPS.jl repo.
    MUMPS.associate_matrix!(mumpsSolver.mumps, H)
    mumps_job!(mumpsSolver.mumps, FACTOR)

    while mumpsSolver.mumps.infog[12] > 0 # while matrix isn't positive definite.
        mumps_job!(mumpsSolver.mumps, FACTOR_CLEANUP)
        τ_old = τ
        τ = max(2 * τ, β)
        for i in axes(H, 1)
            H[i, i] += τ - τ_old # now try H + τ*I
        end
        MUMPS.associate_matrix!(mumpsSolver.mumps, H)

        # TODO better error handling, so the user doesn't have to look up arcane
        # error codes in the MUMPS user manual.
        mumps_job!(mumpsSolver.mumps, FACTOR)
        if mumpsSolver.mumps.infog[12] == 0
            @debug "nonsingular with τ = $τ"
        end
    end
end

function solve!(mumpsSolver::MUMPSHessianSolver, b::Vector{Float64})
    MUMPS.associate_rhs!(mumpsSolver.mumps, b)
    mumps_job!(mumpsSolver.mumps, SOLVE)
    if mumpsSolver.mumps.infog[1] < 0
        throw(ErrorException("MUMPS solve failed with error code "*
                "$(mumpsSolver.mumps.infog[1])"))
    end
    # TODO: is it safe to reuse the same memory here?
    MUMPS.mumps_solve!(b, mumpsSolver.mumps)
end

function cleanup!(mumpsSolver::MUMPSHessianSolver)
    MUMPS.finalize(mumpsSolver.mumps)
    return
end