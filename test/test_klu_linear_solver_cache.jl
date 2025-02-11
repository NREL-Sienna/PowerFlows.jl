using SparseArrays
using LinearAlgebra

function testSolve(k::PF.KLULinSolveCache, A::SparseMatrixCSC{Float64, Int64})
    xb = randn(size(A,1))
    b = deepcopy(xb)
    PF.solve!(k, xb)
    @test A*xb ≈ b
end


@testset "klu linear solver cache" begin
    N = 50
    A = sprandn(Float64, N, N, 0.1)
    while abs(det(A)) < eps(Float64)
        A = sprandn(Float64, N, N, 0.1)
    end
    @assert abs(det(A)) > eps(Float64)

    k = PF.KLULinSolveCache(A)
    PF.full_factor!(k, A)

    testSolve(k, A) # solve is solution.
    
    # solve! is non-allocating.
    b = randn(N)
    temp = @allocated PF.solve!(k, b)
    @test iszero(temp)

    # numeric_refactor! is non allocating, and solves correctly.
    A2 = deepcopy(A)
    copyto!(A2.nzval, rand(length(A2.nzval)))
    temp = @allocated PF.numeric_refactor!(k, A2)
    @test iszero(temp)
    testSolve(k, A2)

    # with reuse_symbolic and check_pattern both true,
    # refactor with a different structure throws.
    B = deepcopy(A)
    n, m = rand(1:N), rand(1:N)
    while (B[n,m] != 0.0)
        n, m = rand(1:N), rand(1:N)
    end
    B[n,m] = 0.1
    @test_throws ArgumentError PF.numeric_refactor!(k, B)
    @test_throws ArgumentError PF.symbolic_refactor!(k, B)

    # symbolic (re)factor and solve. this does allocate.
    PF.symbolic_factor!(k, B)
    PF.numeric_refactor!(k, B)
    testSolve(k, B)

    # automatically refactor if reuse_symbolic is false.
    autoRefactor = PF.KLULinSolveCache(A, false, false)
    @test PF.symbolic_factor!(autoRefactor, B) isa Any # shouldn't throw.
    
    # error handling: singular.
    sing = sparse([1], [1], [0.1], 2, 2)
    # If I add KLU.kluerror(cache.K.common) to the constructor and symbolic_factor!
    # in cache, it still doesn't throw until you call numeric_refactor!. Strange.
    singCache = PF.KLULinSolveCache(sing)
    PF.symbolic_factor!(singCache, sing) 
    @test_throws SingularException PF.numeric_refactor!(singCache, sing)

    # error handling: non-square.
    nonSquare = sprand(10, 11, 0.5)
    @test_throws ArgumentError PF.KLULinSolveCache(nonSquare)

    # error handling: mismtached dimensions.
    b = rand(N-1)
    @test_throws DimensionMismatch PF.solve!(k, b)
end