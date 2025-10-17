struct LCCConverterParameters
    bus::Vector{Int} # bus indices, not bus numbers.
    tap::Matrix{Float64}
    thyristor_angle::Matrix{Float64}
    phi::Matrix{Float64}
    transformer_reactance::Vector{Float64}
    min_thyristor_angle::Vector{Float64}
end

LCCConverterParameters(n_timesteps::Int, n_lccs::Int) = LCCConverterParameters(
    zeros(Int, n_lccs),
    zeros(Float64, n_lccs, n_timesteps),
    zeros(Float64, n_lccs, n_timesteps),
    zeros(Float64, n_lccs, n_timesteps),
    zeros(Float64, n_lccs),
    zeros(Float64, n_lccs),
)

# TODO could add an lcc arc axis/arc lookup pair here and switch to matrices instead of vectors of dicts.
# or better yet, use an ArcAdmittanceMatrix from PNM. But this is fine for now: should only be a few LCCs.
# TODO switch to uniformly using bus numbers. Use lookups where we need indices.
struct LCCParameters
    # bus number pairs, not indices. needed in post-processing for result reporting.
    arcs::Vector{Tuple{Int, Int}}
    # (bus_rectifier, bus_inverter) index pair -> (y_ff, y_tt) admittance pair.
    # single dictionary, reused between time steps.
    branch_admittances::Dict{Tuple{Int, Int}, Tuple{ComplexF64, ComplexF64}}
    # each dict stores (bus_rectifier, bus_inverter) -> power flow.
    # keys here are bus indices, not bus numbers.
    arc_activepower_flow_from_to::Vector{Dict{Tuple{Int, Int}, Float64}}
    arc_activepower_flow_to_from::Vector{Dict{Tuple{Int, Int}, Float64}}
    arc_reactivepower_flow_from_to::Vector{Dict{Tuple{Int, Int}, Float64}}
    arc_reactivepower_flow_to_from::Vector{Dict{Tuple{Int, Int}, Float64}}
    setpoint_at_rectifier::Vector{Bool}
    p_set::Matrix{Float64}
    i_dc::Matrix{Float64}
    dc_line_resistance::Vector{Float64}
    rectifier::LCCConverterParameters
    inverter::LCCConverterParameters
end

LCCParameters(n_timesteps::Int, n_lccs::Int) = LCCParameters(
    sizehint!(Vector{Tuple{Int, Int}}(), n_lccs),
    Dict{Tuple{Int, Int}, Tuple{ComplexF64, ComplexF64}}(),
    [Dict{Tuple{Int, Int}, Float64}() for _ in 1:n_timesteps],
    [Dict{Tuple{Int, Int}, Float64}() for _ in 1:n_timesteps],
    [Dict{Tuple{Int, Int}, Float64}() for _ in 1:n_timesteps],
    [Dict{Tuple{Int, Int}, Float64}() for _ in 1:n_timesteps],
    falses(n_lccs),
    zeros(Float64, n_lccs, n_timesteps),
    zeros(Float64, n_lccs, n_timesteps),
    zeros(Float64, n_lccs),
    LCCConverterParameters(n_timesteps, n_lccs),
    LCCConverterParameters(n_timesteps, n_lccs),
)
