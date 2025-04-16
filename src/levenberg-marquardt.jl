function _newton_powerflow(
    pf::ACPowerFlow{LevenbergMaquardtACPowerFlow},
    data::ACPowerFlowData,
    time_step::Int64;
    kwargs...)
    # begin copy-paste from powerflow_method.jl
    # TODO: find way to reduce code repetition.
    residual = ACPowerFlowResidual(data, time_step)
    x0 = calculate_x0(data, time_step)
    residual(x0, time_step)
    if norm(residual.Rv, Inf) < get(kwargs, :tol, DEFAULT_NR_TOL)
        return true # starting point is already a solution.
    end
    # usually would define J here, but we'll store J inside LevenbergMaquardtData.
    if sum(abs, residual.Rv) > WARN_LARGE_RESIDUAL * length(residual.Rv)
        lg_res, ix = findmax(residual.Rv)
        lg_res_rounded = round(lg_res; sigdigits = 3)
        pow_type = ix % 2 == 1 ? "active" : "reactive"
        bus_ix = div(ix + 1, 2)
        bus_no = axes(data.power_network_matrix, 1)[bus_ix]
        @warn "Initial guess provided results in a large initial residual of $lg_res_rounded. " *
              "Largest residual at bus $bus_no ($bus_ix by matrix indexing; $pow_type power)"
    end
    # end copy-paste from powerflow_method.jl

    lmd = LevenbergMaquardtData(data, time_step)
    linSolveCache = KLULinSolveCache(lmd.A)
    symbolic_factor!(linSolveCache, lmd.A)

    converged, i = _run_powerflow_method(
        time_step,
        x0,
        lmd,
        linSolveCache,
        residual;
        kwargs...,
    )

    if converged
        @info("The LevenbergMaquardtACPowerFlow solver converged after $i iterations.")
        if data.calculate_loss_factors
            residual(x0, time_step)
            lmd.J(time_step)
            calculate_loss_factors(data, lmd.J.Jv, time_step)
        end
        return true
    end
    @error("The LevenbergMaquardtACPowerFlow solver failed to converge.")
    return false
end

"""Implementation of Levenberg-Maquardt."""
struct LevenbergMaquardtData
    J::ACPowerFlowJacobian
    A::SparseMatrixCSC{Float64, Int32} # J' * J + c * I
    b::Vector{Float64}
end

function LevenbergMaquardtData(data::ACPowerFlowData, time_step::Int)
    J = ACPowerFlowJacobian(data, time_step)
    A = _create_JT_J_sparse_structure(data, time_step)
    b = zeros(size(A, 1))
    return LevenbergMaquardtData(J, A, b)
end

"""Run power flow with Levenberg-Maquardt method."""
function _run_powerflow_method(
    time_step::Int,
    x::Vector{Float64},
    lmd::LevenbergMaquardtData,
    linSolveCache::KLULinSolveCache,
    residual::ACPowerFlowResidual;
    kwargs...,
)
    # args compared to other newton-type methods: no stateVector and replace J with lmd.
    maxIterations::Int = get(kwargs, :maxIterations, DEFAULT_NR_MAX_ITER)
    λ::Float64 = get(kwargs, :λ_0, DEFAULT_λ_0)
    tol::Float64 = get(kwargs, :tol, DEFAULT_NR_TOL)
    maxTestλs::Int = get(kwargs, :maxTestλs, DEFAULT_MAX_TEST_λs)
    i, converged = 0, false
    while i < maxIterations && !converged && !isnan(λ)
        λ = update!(lmd, linSolveCache, x, residual, λ, time_step, maxTestλs)
        converged = !isnan(λ) && norm(residual.Rv, Inf) < tol
        i += 1
    end
    return converged, i
end

"""Updates x following to Levenberg-Maquardt: returns true for success,
false for failure. Current procedure for adjusting λ:
(1) improvement with current λ => λ /= DAMPING_DECR and x += Δx.
(2) else, λ *= DAMPING_INCR.
    (2a) improvement with this λ => keep λ the same and x += Δx.
    (2b) else, go back to (2).
Here, "improvement with λ" means: if let `Δx` be solution to
`(J' * J + λ I) Δx = J'*F(x)`, then `norm(F(x+Δx), 2) < norm(F(x), 2)`.
"""
function update!(lmd::LevenbergMaquardtData,
    linSolveCache::KLULinSolveCache,
    x::Vector{Float64},
    residual::ACPowerFlowResidual,
    λ::Float64,
    time_step::Int,
    maxTestλs::Int,
)
    # set lmd.A to J' * J
    lmd.A.nzval .= 0.0
    residual(x, time_step)
    lmd.J(time_step)
    A_plus_eq_BT_B!(lmd.A, lmd.J.Jv)

    # initialize stuff.
    j = 0
    λ_prev = 0.0
    lmd.b .= lmd.J.Jv' * residual.Rv
    @assert !all(lmd.b .== 0.0) "Levenberg-Maquardt is solving `A * Δx = 0`: " *
                                "F(x) is exactly in the kernel of J'(x). Highly degenerate system or a bug."
    # PERF: normalize or re-scale, else the numerically big (power) components
    # will dominate.
    residualSize = LinearAlgebra.dot(residual.Rv, residual.Rv)
    while j < maxTestλs
        # set lmd.A to J' * J + λ * I
        for i in axes(lmd.A, 1)
            lmd.A[i, i] += λ - λ_prev
        end
        numeric_refactor!(linSolveCache, lmd.A)

        Δx = deepcopy(lmd.b)
        solve!(linSolveCache, Δx)
        Δx *= -1
        residual(x + Δx, time_step)
        # PERF: normalize or re-scale
        newResidualSize = LinearAlgebra.dot(residual.Rv, residual.Rv)
        if newResidualSize < residualSize
            x .+= Δx
            if j == 0
                λ /= DAMPING_DECR
            end
            return λ
        end
        λ_prev = λ
        λ *= DAMPING_INCR
        j += 1
    end
    @error "Unable to improve: gave up after increasing damping factor $maxTestλs times."
    return NaN
end

"""Does `A += B' * B`, in a way that preserves the sparse structure of `A`, if possible.
A workaround for the fact that Julia seems to run `dropzeros!(A)` automatically if I just 
do `A .+= B' * B`."""
function A_plus_eq_BT_B!(A::SparseMatrixCSC, B::SparseMatrixCSC)
    M = B' * B
    @assert M.colptr == A.colptr && M.rowval == A.rowval
    A.nzval .+= M.nzval
    return
end

"""Create the sparse structure of `J' * J`. Structurally non zero 2x2 blocks
correspond to ordered pairs of buses with a neighbor in common (or are the same:
i.e. diagonal blocks are all strucutrally nonzero)."""
function _create_JT_J_sparse_structure(data::ACPowerFlowData, time_step::Int)
    # J' * J is dot products of pairs of columns.
    # so look at pairs of columns and check if there's a row in which both are nonzero.
    # i.e. look at pairs of buses and see if they have a neighbor in common.
    rows = Int32[]      # I
    columns = Int32[]   # J
    values = Float64[]  # V

    # an over-estimate: this counts directed paths of 2 edges. Might be many 
    # paths of length 2 between the same pair of points that share a neighbor.
    numEdgePairs = sum(x -> length(x)^2, get_branch_lookup(data); init = 0)
    sizehint!(rows, 4 * numEdgePairs)
    sizehint!(columns, 4 * numEdgePairs)
    sizehint!(values, 4 * numEdgePairs)

    visited = Set{Tuple{Int, Int}}()
    # instead of checking pairs of buses O(n^2 d), we travel 2 hops from each bus O(n d^2)
    # (d = degree of typical vertex, n = number of vertices; d << n.)
    for (ind, nbrs) in enumerate(get_neighbor(data))
        for nbr in nbrs
            for nbrOfNbr in get_neighbor(data)[nbr]
                if !((ind, nbrOfNbr) in visited)
                    bus_from, bus_to = ind, nbrOfNbr
                    for (i, j) in Iterators.product((2 * bus_from - 1, 2 * bus_from),
                        (2 * bus_to - 1, 2 * bus_to))
                        push!(rows, i)
                        push!(columns, j)
                        push!(values, 0.0)
                    end
                    push!(visited, (ind, nbrOfNbr))
                end
            end
        end
    end
    n = size(get_bus_type(data), 1)
    @info "J' * J has sparsity of $(size(values, 1)/(2n)^2)"
    return SparseArrays.sparse(rows, columns, values)
end
