function coordSet(A::SparseMatrixCSC)
    (row, col, _) = SparseArrays.findnz(A)
    return Set(zip(row, col))
end

function testSameStructure(A1::SparseMatrixCSC, A2::SparseMatrixCSC)
    @test A1.rowval == A2.rowval && A1.colptr == A2.colptr
end

@testset "A_plus_eq_BT_B!" begin
    # test basic numeric behavior.
    B = sprand(100, 100, 0.01)
    A = B' * B
    A.nzval .= rand(size(A.nzval, 1))
    A_old = deepcopy(A)
    PF.A_plus_eq_BT_B!(A, B)
    @test isapprox(A_old + B' * B, A)
    testSameStructure(A_old, A)
    # change the sparse structure of A and check it throws an error.
    coordSetA = coordSet(A)
    i, j = rand(1:100), rand(1:100)
    while (i, j) in coordSetA
        i, j = rand(1:100), rand(1:100)
    end
    A_old[i, j] = 1.0
    @test_throws AssertionError PF.A_plus_eq_BT_B!(A_old, B)
    # test that zeros aren't dropped.
    # find a column of B with some nonzero values, and make them all zero.
    firstCol = 1
    while all(A[:, firstCol] .== 0.0)
        firstCol += 1
    end
    vals = SparseArrays.nonzeros(B)
    for i in SparseArrays.nzrange(B, firstCol)
        vals[i] = 0.0
    end
    temp = B' * B
    @assert (firstCol, firstCol) in coordSet(temp)
    @assert temp[firstCol, firstCol] == 0.0
    @assert (firstCol, firstCol) in coordSet(A)
    A[firstCol, firstCol] = 0.0
    @assert (firstCol, firstCol) in coordSet(A)
    A_old = deepcopy(A)
    PF.A_plus_eq_BT_B!(A, B)
    @test (firstCol, firstCol) in coordSet(A)
    testSameStructure(A_old, A)
end

"""Construct a system from a graph: `n` is the number of vertices and 
`edges` the list of edges."""
function createSystemWithTopology(n::Int, edges::Vector{Tuple{Int, Int}})
    @assert n >= 3
    sys = System(100.0)
    for i in 1:n
        if i == 1
            bt = ACBusTypes.REF
        elseif i % 2 == 1
            bt = ACBusTypes.PV
        else
            bt = ACBusTypes.PQ
        end
        b = ACBus(;
            number = i,
            name = "bus$i",
            bustype = bt,
            angle = 0.0,
            magnitude = 1.1,
            voltage_limits = (0.0, 2.0),
            base_voltage = 100,
        )
        add_component!(sys, b)
    end

    for (i, tp) in enumerate(edges)
        (from_ind, to_ind) = tp
        bus_from = get_component(PSY.ACBus, sys, "bus$from_ind")
        bus_to = get_component(PSY.ACBus, sys, "bus$to_ind")
        line = Line(;
            name = "line$i",
            available = true,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = Arc(; from = bus_from, to = bus_to),
            r = 0.00281, # Per-unit
            x = 0.0281, # Per-unit
            b = (from = 0.00356, to = 0.00356), # Per-unit
            rating = 2.0, # Line rating of 200 MVA / System base of 100 MVA
            angle_limits = (min = -0.7, max = 0.7),
        )
        add_component!(sys, line)
    end
    return sys
end

@testset "J' J sparse structure" begin
    # basic types of graphs to test:
    # a chain, a loop, disconnected, dense, and something mildly complex.

    chain_graph = (5, [(1, 2), (2, 3), (3, 4), (4, 5)])
    loop_graph = (5, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])
    denseEdges = Vector{Tuple{Int, Int}}()
    for (i, j) in zip(1:5, 1:5)
        if i < j
            push!(denseEdges, (i, j))
        end
    end
    dense_graph = (5, denseEdges)
    multi_component = (6, [(1, 2), (2, 3), (4, 5)])
    # start with a loop of size 9, add in a few other edges
    complex_graph_edges = Vector{Tuple{Int, Int}}([(i, i + 1) for i in 1:8])
    push!(complex_graph_edges, (9, 1))
    append!(complex_graph_edges, [(3, 9), (1, 4), (3, 5), (6, 8)])
    complex_graph = (9, complex_graph_edges)

    graphs = [chain_graph, dense_graph, multi_component, complex_graph]
    for graph in graphs
        sys = createSystemWithTopology(graph...)
        pf = ACPowerFlow()
        data = PF.PowerFlowData(pf, sys; check_connectivity = false)
        # no need to actually initialize J.Jv: we only care about the structure.
        J = PF.ACPowerFlowJacobian(data, 1)
        JTJ_structure = PF._create_JT_J_sparse_structure(data, 1)
        M = J.Jv' * J.Jv
        testSameStructure(M, JTJ_structure)
    end
end
