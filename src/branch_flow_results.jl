
const BranchFlowEntry = @NamedTuple{
    name::String,
    bus_from::Int,
    bus_to::Int,
    P_from_to::Float64,
    P_to_from::Float64,
    P_losses::Float64,
    Q_from_to::Float64,
    Q_to_from::Float64,
    Q_losses::Float64,
}

const _COL_P_FROM_TO = 1
const _COL_P_TO_FROM = 2
const _COL_P_LOSSES = 3
const _COL_Q_FROM_TO = 4
const _COL_Q_TO_FROM = 5
const _COL_Q_LOSSES = 6
const _N_FLOW_COLS = 6

"""Column-major storage for branch flow results. Entries are added via `push!` using
`BranchFlowEntry` named tuples for readability; flow values are stored in a pre-allocated
matrix whose columns can be extracted directly for DataFrame construction."""
mutable struct BranchFlowResults
    names::Vector{String}
    bus_from::Vector{Int}
    bus_to::Vector{Int}
    flows::Matrix{Float64}   # n × _N_FLOW_COLS
    count::Int
end

function BranchFlowResults(n::Int)
    return BranchFlowResults(
        Vector{String}(undef, n),
        Vector{Int}(undef, n),
        Vector{Int}(undef, n),
        Matrix{Float64}(undef, n, _N_FLOW_COLS),
        0,
    )
end

function Base.push!(r::BranchFlowResults, e::BranchFlowEntry)
    i = r.count += 1
    r.names[i] = e.name
    r.bus_from[i] = e.bus_from
    r.bus_to[i] = e.bus_to
    r.flows[i, _COL_P_FROM_TO] = e.P_from_to
    r.flows[i, _COL_P_TO_FROM] = e.P_to_from
    r.flows[i, _COL_P_LOSSES] = e.P_losses
    r.flows[i, _COL_Q_FROM_TO] = e.Q_from_to
    r.flows[i, _COL_Q_TO_FROM] = e.Q_to_from
    r.flows[i, _COL_Q_LOSSES] = e.Q_losses
    return r
end
