"""In order to in-place modify the numeric values of a CHOLMOD matrix, we need to 
write our own wrapper around `CHOLMOD.Sparse`. """
struct FixedStructureCHOLMOD{T, I} <: SparseArrays.AbstractSparseMatrix{T, I}
    _mat::SparseArrays.CHOLMOD.Sparse{T, I}
    _values::Vector{T}
end

Base.show(mat::FixedStructureCHOLMOD) = print("FixedStructureCHOLMOD with $(size(mat)) and $(length(mat._values)) values")

function FixedStructureCHOLMOD(A::SparseMatrixCSC{T, I}) where {T<:VTypes, I<:ITypes}
    cholmod_mat = SparseArrays.CHOLMOD.Sparse(LinearAlgebra.Symmetric(A))

    chol = unsafe_load(cholmod_mat.ptr)
    x_ptr = chol.x
    # max capacity may not be the same as the number of nonzeros:
    # instead of chol.nzmax, we want chol.p[chol.ncol]
    # (what about packed matrices, though?)
    nonzero_count = unsafe_load(Ptr{I}(chol.p), chol.ncol + 1)
    values = unsafe_wrap(Vector{T}, Ptr{T}(x_ptr), nonzero_count)
    return FixedStructureCHOLMOD(cholmod_mat, values)
end

"""
    set_values!(mat::FixedStructureCHOLMOD, new_vals::AbstractVector{Float64})
In-place update of the numeric values in the CHOLMOD matrix.
"""
function set_values!(mat::FixedStructureCHOLMOD{T, I}, new_vals::AbstractVector{T}) where {T<:VTypes, I<:ITypes}
    @assert length(new_vals) == length(mat._values) "New values must be same length as allocated storage!"
    copyto!(mat._values, new_vals)
    return mat
end

function symbolic_factor(mat::FixedStructureCHOLMOD{T, I}) where {T<:VTypes, I<:ITypes}
    # I'm really feeding in a CHOLMOD_PATTERN object here, ie xtype = 0, but Julia's
    # interface doesn't seem to support that.
    return SparseArrays.CHOLMOD.symbolic(mat._mat)
end

function numeric_factor!(F::SparseArrays.CHOLMOD.Factor{T, I}, mat::FixedStructureCHOLMOD{T, I}) where {T<:VTypes, I<:ITypes}
    SparseArrays.CHOLMOD.factorize!(mat._mat, F)
    return
end

Base.size(mat::FixedStructureCHOLMOD) = size(mat._mat)
Base.getindex(mat::FixedStructureCHOLMOD, I...) = getindex(mat._mat, I...)
Base.setindex!(mat::FixedStructureCHOLMOD, v, I...) = setindex!(mat._mat, v, I...)
Base.eltype(mat::FixedStructureCHOLMOD) = eltype(mat._mat)
