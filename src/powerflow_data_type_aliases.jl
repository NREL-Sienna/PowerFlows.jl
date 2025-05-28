"""A type alias for a [`PowerFlowData`](@ref) struct whose type parameters
are configured for the [`ACPowerFlow`](@ref) method."""
const ACPowerFlowData = PowerFlowData{
    PNM.Ybus{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    },
    Nothing,
}

"""A type alias for a [`PowerFlowData`](@ref) struct whose type parameters
are configured for the [`PTDFDCPowerFlow`](@ref) method ."""
const PTDFPowerFlowData = PowerFlowData{
    PNM.PTDF{
        Tuple{Vector{Int64}, Vector{String}},
        Tuple{Dict{Int64, Int64}, Dict{String, Int64}},
        Matrix{Float64},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

"""A type alias for a [`PowerFlowData`](@ref)  struct whose type parameters
are configured for the [`vPTDFDCPowerFlow`](@ref) method."""
const vPTDFPowerFlowData = PowerFlowData{
    PNM.VirtualPTDF{
        Tuple{Vector{String}, Vector{Int64}},
        Tuple{Dict{String, Int64}, Dict{Int64, Int64}},
    },
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
}

"""A type alias for a [`PowerFlowData`](@ref)  struct whose type parameters
are configured for the [`DCPowerFlow`](@ref) method."""
const ABAPowerFlowData = PowerFlowData{
    PNM.ABA_Matrix{
        Tuple{Vector{Int64}, Vector{Int64}},
        Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
        PNM.KLU.KLUFactorization{Float64, Int64},
    },
    PNM.BA_Matrix{
        Tuple{Vector{Int64}, Vector{String}},
        Tuple{Dict{Int64, Int64}, Dict{String, Int64}}},
}
