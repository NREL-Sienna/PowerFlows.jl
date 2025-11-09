module PowerFlowsCUDAExt

using PowerFlows
using CUDA
using CUDA: CUSOLVER, CUSPARSE
using SparseArrays
using LinearAlgebra

import PowerFlows: LinearSolverCache, symbolic_factor!, symbolic_refactor!,
                   numeric_refactor!, solve!, solve_w_refinement,
                   get_reuse_symbolic, full_factor!,
                   DEFAULT_REFINEMENT_MAX_ITER

"""
CUDA-based linear solver using cusolverSpDcsrlsvlu for sparse linear systems.
This provides a GPU-accelerated alternative to KLU for Newton-Raphson power flow.

Note: This extension is automatically loaded when CUDA.jl is available.
"""

"""A cached linear solver using cuSOLVER's cusolverSpDcsrlsvlu.
Uses GPU acceleration for sparse LU factorization and solve.

# Fields:
- `n::Int`: matrix dimension
- `nnz::Int`: number of non-zeros
- `handle`: cuSOLVER handle
- `d_A`: device storage for matrix values
- `d_csrRowPtr`: device storage for CSR row pointers
- `d_csrColInd`: device storage for CSR column indices
- `d_x`: device storage for solution vector
- `d_b`: device storage for right-hand side
- `h_csrRowPtr`: host storage for row pointers (cached)
- `h_csrColInd`: host storage for column indices (cached)
- `reuse_symbolic::Bool`: reuse the symbolic factorization
- `check_pattern::Bool`: verify sparsity pattern matches on refactor
- `pivot_threshold::Float64`: pivoting threshold for stability
- `info`: cuSOLVER info structure
"""
mutable struct CUSOLVERLinSolveCache{T<:Integer} <: PowerFlows.LinearSolverCache{T}
    n::Int
    nnz::Int
    handle::Any  # Will be CUSOLVER.cusolverSpHandle
    d_A::Any     # Will be CuVector{Float64}
    d_csrRowPtr::Any  # Will be CuVector{T}
    d_csrColInd::Any  # Will be CuVector{T}
    d_x::Any     # Will be CuVector{Float64}
    d_b::Any     # Will be CuVector{Float64}
    h_csrRowPtr::Vector{T}
    h_csrColInd::Vector{T}
    reuse_symbolic::Bool
    check_pattern::Bool
    pivot_threshold::Float64
    info::Any    # Will be cusolverSpInfo

    function CUSOLVERLinSolveCache{T}(
        n::Int,
        nnz::Int,
        handle,
        d_A,
        d_csrRowPtr,
        d_csrColInd,
        d_x,
        d_b,
        h_csrRowPtr::Vector{T},
        h_csrColInd::Vector{T},
        reuse_symbolic::Bool,
        check_pattern::Bool,
        pivot_threshold::Float64,
        info,
    ) where {T<:Integer}
        new{T}(n, nnz, handle, d_A, d_csrRowPtr, d_csrColInd, d_x, d_b,
               h_csrRowPtr, h_csrColInd, reuse_symbolic, check_pattern,
               pivot_threshold, info)
    end
end

"""Constructor for CUSOLVERLinSolveCache.
Converts CSC matrix to CSR format and allocates GPU memory.
"""
function CUSOLVERLinSolveCache(
    A::SparseMatrixCSC{Float64, T};
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true,
    pivot_threshold::Float64 = 1.0,
) where {T<:Integer}
    n, m = size(A)
    if n != m
        throw(ArgumentError("Matrix must be square. Got size ($n, $m)."))
    end

    # Convert CSC to CSR (transpose to get CSR from CSC)
    A_csr = SparseMatrixCSC(A')
    nnz = length(A_csr.nzval)

    # Create cuSOLVER handle
    handle = CUSOLVER.cusolverSpCreate()

    # Allocate device memory
    d_A = CUDA.CuVector{Float64}(undef, nnz)
    d_csrRowPtr = CUDA.CuVector{T}(undef, n + 1)
    d_csrColInd = CUDA.CuVector{T}(undef, nnz)
    d_x = CUDA.CuVector{Float64}(undef, n)
    d_b = CUDA.CuVector{Float64}(undef, n)

    # Store host copies for pattern checking
    h_csrRowPtr = Vector{T}(A_csr.colptr)
    h_csrColInd = Vector{T}(A_csr.rowval)

    # Create info structure
    info = CUSOLVER.cusolverSpCreateCsrluInfoHost()

    cache = CUSOLVERLinSolveCache{T}(
        n, nnz, handle, d_A, d_csrRowPtr, d_csrColInd, d_x, d_b,
        h_csrRowPtr, h_csrColInd, reuse_symbolic, check_pattern,
        pivot_threshold, info
    )

    return cache
end

"""Get whether symbolic factorization is reused."""
get_reuse_symbolic(cache::CUSOLVERLinSolveCache) = cache.reuse_symbolic

"""Performs symbolic factorization of matrix A and stores it in cache."""
function symbolic_factor!(
    cache::CUSOLVERLinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T<:Integer}
    if !(size(A, 1) == cache.n && size(A, 2) == cache.n)
        throw(DimensionMismatch(
            "Can't factor: matrix has different dimensions."
        ))
    end

    # Convert CSC to CSR
    A_csr = SparseMatrixCSC(A')

    # Update cached pattern
    cache.h_csrRowPtr = Vector{T}(A_csr.colptr)
    cache.h_csrColInd = Vector{T}(A_csr.rowval)

    # Copy to device
    CUDA.copyto!(cache.d_csrRowPtr, cache.h_csrRowPtr)
    CUDA.copyto!(cache.d_csrColInd, cache.h_csrColInd)

    # For cusolverSpDcsrlsvlu, symbolic factorization is done implicitly
    # during the first numeric factorization, so we just prepare the structure

    return nothing
end

"""Symbolic refactorization with optional pattern checking."""
function symbolic_refactor!(
    cache::CUSOLVERLinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T<:Integer}
    if cache.reuse_symbolic && cache.check_pattern
        if !(size(A, 1) == cache.n && size(A, 2) == cache.n)
            throw(DimensionMismatch(
                "Can't refactor: new matrix has different dimensions."
            ))
        end

        # Check pattern matches
        A_csr = SparseMatrixCSC(A')
        new_rowptr = Vector{T}(A_csr.colptr)
        new_colind = Vector{T}(A_csr.rowval)

        if new_rowptr != cache.h_csrRowPtr || new_colind != cache.h_csrColInd
            throw(ArgumentError(
                "Matrix has different sparse structure. Either create cache with " *
                "reuse_symbolic = false, or call symbolic_factor! instead."
            ))
        end
    elseif !cache.reuse_symbolic
        symbolic_factor!(cache, A)
    end

    return nothing
end

"""Performs numeric factorization of matrix A."""
function numeric_refactor!(
    cache::CUSOLVERLinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T<:Integer}
    if cache.check_pattern
        A_csr = SparseMatrixCSC(A')
        new_rowptr = Vector{T}(A_csr.colptr)
        new_colind = Vector{T}(A_csr.rowval)

        if new_rowptr != cache.h_csrRowPtr || new_colind != cache.h_csrColInd
            throw(ArgumentError(
                "Cannot numeric_refactor: matrix has different sparse structure."
            ))
        end
    end

    # Convert to CSR and copy values to device
    A_csr = SparseMatrixCSC(A')
    CUDA.copyto!(cache.d_A, A_csr.nzval)

    # cusolverSpDcsrlsvlu performs factorization during solve
    # For now, we just update the matrix values

    return nothing
end

"""Solves the linear system using cusolverSpDcsrlsvlu.
Modifies B in-place with the solution."""
function solve!(
    cache::CUSOLVERLinSolveCache{T},
    B::StridedVecOrMat{Float64},
) where {T<:Integer}
    if size(B, 1) != cache.n
        throw(DimensionMismatch(
            "Need size(B, 1) to equal $(cache.n), but got $(size(B, 1))."
        ))
    end

    if stride(B, 1) != 1
        throw(ArgumentError("B must have unit strides"))
    end

    # Handle vector case
    if B isa Vector
        # Copy RHS to device
        CUDA.copyto!(cache.d_b, B)

        # Create matrix descriptor for CSR format
        # Note: cuSOLVER expects 0-based indexing, so we need to adjust
        descrA = CUSPARSE.CuSparseMatrixDescriptor()
        CUSPARSE.cusparseSetMatIndexBase(descrA, CUSPARSE.CUSPARSE_INDEX_BASE_ZERO)
        CUSPARSE.cusparseSetMatType(descrA, CUSPARSE.CUSPARSE_MATRIX_TYPE_GENERAL)

        # Adjust indices for 0-based indexing
        d_csrRowPtr_zero = cache.d_csrRowPtr .- T(1)
        d_csrColInd_zero = cache.d_csrColInd .- T(1)

        # Solve using cusolverSpDcsrlsvlu
        tol = eps(Float64)
        reorder = 1  # Enable reordering for better stability
        singularity = Ref{Cint}(0)

        CUSOLVER.cusolverSpDcsrlsvluHost(
            cache.handle,
            cache.n,
            cache.nnz,
            descrA,
            cache.d_A,
            d_csrRowPtr_zero,
            d_csrColInd_zero,
            cache.d_b,
            tol,
            reorder,
            cache.d_x,
            singularity
        )

        if singularity[] >= 0
            @warn "Matrix appears to be singular at position $(singularity[])"
        end

        # Copy solution back to host
        CUDA.copyto!(B, cache.d_x)

    else
        # Handle matrix case (multiple RHS)
        for i in 1:size(B, 2)
            b_col = view(B, :, i)
            solve!(cache, b_col)
        end
    end

    return B
end

"""Solve with iterative refinement for better accuracy on ill-conditioned systems."""
function solve_w_refinement(
    cache::CUSOLVERLinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
    B::StridedVecOrMat{Float64},
    tol::Float64 = 1e-6,
) where {T<:Integer}
    bNorm = norm(B, 1)
    XB = zeros(size(B))
    r = B - A * XB
    MAX_ITERS = DEFAULT_REFINEMENT_MAX_ITER
    iters = 0

    while iters < MAX_ITERS && norm(r, 1) >= bNorm * tol
        lastError = norm(r, 1)
        solve!(cache, r)
        XB .+= r
        r .= B - A * XB
        iters += 1

        if norm(r, 1) > lastError
            @error "Iterative refinement failed: error is getting worse."
            return XB
        end
    end

    @debug "Iterative refinement converged in $iters iterations."
    return XB
end

"""Cleanup function to free cuSOLVER resources."""
function Base.finalize(cache::CUSOLVERLinSolveCache)
    try
        CUSOLVER.cusolverSpDestroy(cache.handle)
        CUSOLVER.cusolverSpDestroyCsrluInfoHost(cache.info)
    catch e
        @warn "Error during CUSOLVERLinSolveCache cleanup" exception = e
    end
end

# Export the CUDA solver type so it can be used
export CUSOLVERLinSolveCache

end # module PowerFlowsCUDAExt
