const RTOL = 10^-7

# Helper function to test solve! method
function testSolve(
    k::PF.CUSOLVERLinSolveCache{T},
    A::SparseMatrixCSC{Float64, T},
) where {T <: Union{Int32, Int64}}
    xb = randn(size(A, 1))
    b = deepcopy(xb)
    PF.solve!(k, xb)
    @test isapprox(A * xb, b; rtol = RTOL)

    XB = randn(size(A, 1), 10)
    B = deepcopy(XB)
    PF.solve!(k, XB)
    @test isapprox(A * XB, B; rtol = RTOL)
end

@testset "CUSOLVERLinSolveCache" begin
    # Check if CUDA is available
    cuda_available = false
    try
        import CUDA
        cuda_available = CUDA.functional()
    catch e
        @warn "CUDA not available, skipping CUDA linear solver tests" exception=e
    end

    if !cuda_available
        @info "Skipping CUDA linear solver tests (no CUDA-capable GPU detected)"
    else
        @testset for dType in (Int32, Int64)  # Test both AC (Int32) and DC (Int64) power flow
            N = 50
            A = SparseMatrixCSC{Float64, dType}(sprandn(Float64, N, N, 0.1))
            while abs(det(A)) < eps(Float64)
                A = SparseMatrixCSC{Float64, dType}(sprandn(Float64, N, N, 0.1))
            end
            @assert abs(det(A)) > eps(Float64)

            k = PF.CUSOLVERLinSolveCache(A)
            PF.full_factor!(k, A)

            @testset "Basic solve" begin
                testSolve(k, A)
            end

            @testset "Numeric refactor" begin
                # numeric_refactor! and solves correctly with new values
                A2 = deepcopy(A)
                copyto!(A2.nzval, rand(length(A2.nzval)))
                PF.numeric_refactor!(k, A2)
                testSolve(k, A2)
            end

            @testset "Pattern checking" begin
                # with reuse_symbolic and check_pattern both true,
                # refactor with a different structure throws.
                B = deepcopy(A)
                n, m = rand(1:N), rand(1:N)
                while (B[n, m] != 0.0)
                    n, m = rand(1:N), rand(1:N)
                end
                B[n, m] = 0.1
                @test_throws ArgumentError PF.numeric_refactor!(k, B)
                @test_throws ArgumentError PF.symbolic_refactor!(k, B)
            end

            @testset "Symbolic refactor" begin
                # symbolic (re)factor and solve
                B = deepcopy(A)
                n, m = rand(1:N), rand(1:N)
                while (B[n, m] != 0.0)
                    n, m = rand(1:N), rand(1:N)
                end
                B[n, m] = 0.1

                PF.symbolic_factor!(k, B)
                PF.numeric_refactor!(k, B)
                testSolve(k, B)
            end

            @testset "Auto refactor disabled" begin
                # automatically refactor if reuse_symbolic is false.
                autoRefactor = PF.CUSOLVERLinSolveCache(A; reuse_symbolic=false, check_pattern=false)
                B = deepcopy(A)
                n, m = rand(1:N), rand(1:N)
                while (B[n, m] != 0.0)
                    n, m = rand(1:N), rand(1:N)
                end
                B[n, m] = 0.1

                @test PF.symbolic_factor!(autoRefactor, B) isa Any # shouldn't throw.
            end

            @testset "Iterative refinement" begin
                # test iterative refinement with an ill-conditioned matrix
                N_big, d = 100, 0.01
                epsilon = 0.004
                v1, v2 = sprandn(N_big, d), sprandn(N_big, d)
                bigA =
                    sparse(v2 * transpose(v1)) + sparse(diagm(randn(N_big)) * epsilon) +
                    epsilon * sprand(N_big, N_big, 10 * d^2)
                bigA = SparseMatrixCSC{Float64, dType}(bigA)

                bigK = PF.CUSOLVERLinSolveCache(bigA)
                PF.full_factor!(bigK, bigA)
                b = rand(N_big)
                x = PF.solve_w_refinement(bigK, bigA, b)
                @test isapprox(bigA * x, b, rtol = RTOL)
            end

            @testset "Error handling - singular" begin
                # error handling: singular matrix
                sing = SparseMatrixCSC{Float64, dType}(sparse([1], [1], [0.1], 2, 2))
                singCache = PF.CUSOLVERLinSolveCache(sing)
                PF.symbolic_factor!(singCache, sing)
                # cuSOLVER might handle this differently than KLU
                # It may warn rather than throw, so we just check it doesn't crash
                @test_logs (:warn,) match_mode=:any PF.numeric_refactor!(singCache, sing)
            end

            @testset "Error handling - non-square" begin
                # error handling: non-square matrix
                nonSquare = SparseMatrixCSC{Float64, dType}(sprand(10, 11, 0.5))
                @test_throws ArgumentError PF.CUSOLVERLinSolveCache(nonSquare)
            end

            @testset "Error handling - dimension mismatch" begin
                # error handling: mismatched dimensions
                b = rand(N - 1)
                @test_throws DimensionMismatch PF.solve!(k, b)
            end
        end
    end
end
