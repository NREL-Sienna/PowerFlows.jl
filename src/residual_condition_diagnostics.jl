#=
Per-iteration solver diagnostics (`log_solver_diagnostics`) and a fold /
voltage-collapse bail-out (`stop_at_fold`).

λ_min is taken on the bus-voltage Schur complement S = A − B·D⁻¹·C of the blocked
Jacobian J = [A B; C D], whose non-bus tail (LCC + VSC + area interchange,
`state_tail_length(data, dcn)` rows/cols) is the trailing block. The (1,1) block of
J⁻¹ is exactly S⁻¹, so v ↦ (J⁻¹·[v; 0])[1:nb] applies S⁻¹ from the *existing*
factorization of J — no second matrix or factorization. With no tail, S = J. The
monitor line and the bail-out share one refactor and one eigensolve via
`run_solver_diagnostics!`.
=#

"""Round to 4 significant figures."""
_sf4(x) = round(x; sigdigits = 4)

"""Applies S⁻¹ via a back-solve of the full `J`: pads `v` with zeros in the
LCC-tail slots, applies `J⁻¹`, returns the leading `n_bus` block."""
struct SchurInverseOperator{C}
    cache::C
    n_bus::Int
    buffer::Vector{Float64}   # padded RHS, length = full state
end

"""Condition estimate κ̂(J), or `NaN` when the backend exposes none. The NaN
fallback is restricted to the non-KLU `PFLinearSolverCache` members so the concrete
`KLULinSolveCache` doesn't shadow the KLU method onto the NaN path."""
_diag_condest(cache::PNM.KLULinSolveCache) = condest!(cache)
_diag_condest(::Union{PNM.AAFactorCache, PardisoLinSolveCache}) = NaN

function (op::SchurInverseOperator)(v::AbstractVector{Float64})
    b = op.buffer
    @inbounds begin
        copyto!(view(b, 1:(op.n_bus)), v)
        fill!(view(b, (op.n_bus + 1):length(b)), 0.0)
    end
    solve!(op.cache, b)
    # KrylovKit stores each returned vector, so hand back a fresh copy of the
    # bus block rather than the reused buffer.
    return b[1:(op.n_bus)]
end

"""Smallest-magnitude eigenvalue of the Schur complement `S` by inverse iteration:
KrylovKit finds the largest-magnitude eigenvalue `μ` of `S⁻¹` and returns `1/μ`.
`S` is non-symmetric, so the result may be complex. Returns `(λ_min, converged)`,
with `converged = false` (and `λ_min = NaN ± NaN im`) on any failure."""
function _schur_min_eigenvalue(
    op::SchurInverseOperator;
    tol::Float64 = 1e-6,
    maxiter::Int = 200,
    krylovdim::Int = 30,
)::Tuple{ComplexF64, Bool}
    n = op.n_bus
    v0 = fill(1.0 / sqrt(n), n)   # deterministic init for reproducible logs
    vals, _, info = KrylovKit.eigsolve(op, v0, 1, :LM; tol, maxiter, krylovdim)
    if info.converged < 1 || isempty(vals)
        return complex(NaN, NaN), false
    end
    μ = vals[1]
    # A (near-)zero dominant eigenvalue of S⁻¹ makes 1/μ overflow; treat it as
    # not-converged rather than reporting an Inf eigenvalue of S.
    if abs(μ) <= eps(Float64)
        return complex(NaN, NaN), false
    end
    return inv(ComplexF64(μ)), true
end

"""Format a (possibly complex) eigenvalue to 4 significant figures as `a` or
`a ± b im`."""
function _fmt_eig(z::Number)
    iz = imag(z)
    return if iz == 0
        "$(_sf4(real(z)))"
    else
        "$(_sf4(real(z))) $(iz < 0 ? "-" : "+") $(_sf4(abs(iz)))im"
    end
end

"""The system bus number for the `bus_ix`-th bus (reduced ordering)."""
_diag_bus_number(data::ACPowerFlowData, bus_ix::Int) =
    axes(data.power_network_matrix, 1)[bus_ix]

const _LCC_RESIDUAL_ROW_NAMES =
    ("P-setpoint", "DC-line balance", "rectifier α-limit", "inverter α-limit")

"""Describe a residual entry that falls in the LCC tail (4 rows per LCC)."""
function _describe_lcc_residual_entry(data::ACPowerFlowData, tail_ix::Int)
    i = div(tail_ix - 1, 4) + 1
    row = mod1(tail_ix, 4)
    from_no, to_no = data.lcc.arcs[i]
    return "LCC $(from_no)→$(to_no) ($(_LCC_RESIDUAL_ROW_NAMES[row]))"
end

"""Describe a residual entry that falls in the area-interchange tail (1 row per
controlled area, keyed by `tail_ix` rather than vector position so a mid-solve
de-enrollment renumbering can't desync this from `_set_area_tail_residuals!`)."""
function _describe_area_residual_entry(data::ACPowerFlowData, tail_ix::Int)
    for area in data.area_interchange.areas
        area.tail_ix == tail_ix && return "area $(area.name) (NI−PDES)"
    end
    error(
        "area_interchange tail_ix=$tail_ix not found among " *
        "$(length(data.area_interchange.areas)) controlled areas",
    )
end

"""`(bus index, 1-based row within that bus's block)` for variable-block
formulations, from the `bus_state_offset` table."""
function _locate_variable_block(offsets::AbstractVector, ix::Int)
    b = searchsortedlast(offsets, ix)
    return b, ix - Int(offsets[b]) + 1
end

# Formulation-aware label for the entry where ‖F‖∞ is attained. The bus block is
# laid out first; the polar-only tail is `[LCC][VSC][area]` (area rows LAST, see
# `area_tail_offset`) — rectangular/mixed never carry an area tail (area-interchange
# control is rejected at construction for those formulations).
function _describe_residual_entry(
    ::ACPowerFlowResidual,
    data::ACPowerFlowData,
    ::Int,
    ix::Int,
)
    n_bus_eqs = 2 * size(data.bus_type, 1)
    if ix <= n_bus_eqs
        bus_ix = div(ix - 1, 2) + 1
        return "bus $(_diag_bus_number(data, bus_ix)) ($(isodd(ix) ? "P" : "Q"))"
    end
    if n_controlled_areas(data) > 0
        area_off = area_tail_offset(data, get_dc_network(data))
        if ix > area_off
            return _describe_area_residual_entry(data, ix - area_off)
        end
    end
    return _describe_lcc_residual_entry(data, ix - n_bus_eqs)
end

function _describe_residual_entry(
    r::ACRectangularCIResidual,
    data::ACPowerFlowData,
    ::Int,
    ix::Int,
)
    if ix <= r.total_bus_state
        b, row = _locate_variable_block(r.bus_state_offset, ix)
        labels = ("ΔI_re", "ΔI_im", "|V|²−V_set²")   # PV uses the 3rd row
        return "bus $(_diag_bus_number(data, b)) ($(labels[row]))"
    end
    return _describe_lcc_residual_entry(data, ix - r.total_bus_state)
end

function _describe_residual_entry(
    r::ACMixedCPBResidual,
    data::ACPowerFlowData,
    time_step::Int,
    ix::Int,
)
    if ix <= r.total_bus_state
        b, row = _locate_variable_block(r.bus_state_offset, ix)
        bt = data.bus_type[b, time_step]
        labels = if bt == PSY.ACBusTypes.PV
            ("ΔP", "|V|²−V_set²")
        elseif bt == PSY.ACBusTypes.PQ
            ("ΔI_im", "ΔI_re")
        else  # REF
            ("ΔI_re", "ΔI_im")
        end
        return "bus $(_diag_bus_number(data, b)) ($(labels[row]))"
    end
    return _describe_lcc_residual_entry(data, ix - r.total_bus_state)
end

# ---------------------------------------------------------------------------
# Fold / voltage-collapse bail-out state and the shared per-iteration hook.
# ---------------------------------------------------------------------------

# UNSEEN = pre-first-observation. A non-finite/non-converged λ_min is deliberately
# not a sign here; the bail decision treats it as a conservative abort.
IS.@scoped_enum(EigvalSign, UNSEEN = 0, NEGATIVE = -1, POSITIVE = 1)

"""Per-solve scratch for [`run_solver_diagnostics!`](@ref): previous ‖F‖∞ (`prev_F`),
last-seen sign of `real(λ_min)` (`eig_sign`), and a reusable padded RHS (`buffer`) so
the Schur operator allocates nothing per iteration."""
mutable struct SolverDiagnosticsState
    prev_F::Float64
    eig_sign::EigvalSign
    buffer::Vector{Float64}
end

SolverDiagnosticsState(n_state::Int) =
    SolverDiagnosticsState(NaN, EigvalSign.UNSEEN, Vector{Float64}(undef, n_state))

"""Set up a solver loop's diagnostics: returns `(monitor, diag_state)`, allocating
the scratch only when a diagnostic or the bail-out is on so the default solve path
allocates nothing. `diag_state` is `nothing` when neither is requested."""
function setup_solver_diagnostics(
    J::Union{ACPowerFlowJacobian, ACRectangularCIJacobian, ACMixedCPBJacobian},
    bail::Bool,
)
    monitor = get_log_solver_diagnostics(J.data)
    diag_state = (monitor || bail) ? SolverDiagnosticsState(size(J.Jv, 1)) : nothing
    return monitor, diag_state
end

"""Update `state.eig_sign` from `λ_min` and decide the fold bail-out. A non-converged
or non-finite `real(λ_min)` is a conservative bail (warn + abort), never a silent
no-op. An exact-zero real part keeps the prior sign. Returns `true` to abort."""
function _decide_eig_sign_switch!(
    state::SolverDiagnosticsState,
    label::AbstractString,
    λ_min::ComplexF64,
    converged::Bool,
)::Bool
    s = real(λ_min)
    if !converged || !isfinite(s)
        @warn "$label: λ_min(S) is indeterminate " *
              "(converged = $converged, λ_min = $(_fmt_eig(λ_min))); treating it as " *
              "a fold / voltage-collapse signature and aborting the search."
        return true
    end
    prev = state.eig_sign
    current = s > 0 ? EigvalSign.POSITIVE : (s < 0 ? EigvalSign.NEGATIVE : prev)
    state.eig_sign = current

    # A switch needs a real, previously-seen sign that differs from the new one.
    switched = prev != EigvalSign.UNSEEN && current != prev
    switched || return false

    @warn "$label: λ_min(S) real part switched sign " *
          "($(prev == EigvalSign.POSITIVE ? "+" : "−") → " *
          "$(current == EigvalSign.POSITIVE ? "+" : "−")), " *
          "λ_min = $(_fmt_eig(λ_min)). This is a fold / voltage-collapse " *
          "signature; aborting the search."
    return true
end

"""Applies `(JᵀJ)⁻¹ = J⁻¹·J⁻ᵀ` via two back-solves against a KLU factorization of
`J` (`J` real ⇒ adjoint = transpose). KLU-only by construction: it needs the transposed
solve `tsolve!` from an existing factorization of `J`, which AppleAccelerate doesn't expose."""
struct GramInverseOperator{C}
    cache::C   # KLU factorization of J
    buffer::Vector{Float64}
end

function (op::GramInverseOperator)(v::AbstractVector{Float64})
    y = copyto!(op.buffer, v)
    tsolve!(op.cache, y)          # Jᵀ y = v  ⇒ y = J⁻ᵀ v
    solve!(op.cache, y)           # J z = y   ⇒ z = J⁻¹ y = (JᵀJ)⁻¹ v
    return copy(y)
end

"""Bus label for a *state*-vector index (mirrors the residual block layout, but
the per-row variable differs, so only the bus is reported)."""
# LCC *state* variables per converter (cf. state_indexing_helpers.jl), as opposed
# to the residual rows (`_LCC_RESIDUAL_ROW_NAMES`). i_dc is precomputed, not a state.
const _LCC_STATE_NAMES = ("tap_r", "tap_i", "α_r", "α_i")

function _describe_lcc_state_entry(data::ACPowerFlowData, tail_ix::Int)
    i = div(tail_ix - 1, 4) + 1
    from_no, to_no = data.lcc.arcs[i]
    return "LCC $(from_no)→$(to_no) ($(_LCC_STATE_NAMES[mod1(tail_ix, 4)]))"
end

function _state_bus_label(::ACPowerFlowResidual, data::ACPowerFlowData, ix::Int)
    n_bus_eqs = 2 * size(data.bus_type, 1)
    ix <= n_bus_eqs && return "bus $(_diag_bus_number(data, div(ix - 1, 2) + 1))"
    return _describe_lcc_state_entry(data, ix - n_bus_eqs)
end

function _state_bus_label(
    r::Union{ACRectangularCIResidual, ACMixedCPBResidual},
    data::ACPowerFlowData,
    ix::Int,
)
    if ix <= r.total_bus_state
        b, _ = _locate_variable_block(r.bus_state_offset, ix)
        return "bus $(_diag_bus_number(data, b))"
    end
    return _describe_lcc_state_entry(data, ix - r.total_bus_state)
end

"""Top `k` entries of `v` by magnitude, as `(label(ix), v[ix])` pairs."""
function _top_entries(v::AbstractVector{Float64}, k::Int, label)
    idx = partialsortperm(v, 1:min(k, length(v)); by = abs, rev = true)
    return [(label(i), round(v[i]; sigdigits = 4)) for i in idx]
end

"""Smallest singular triplet `(σ_min, v_min, u_min)` of the AC power flow Jacobian
`J` at `x0` (default: flat start), plus the `residual`/`jac` objects so callers can
reuse the labeling helpers. `v_min` is the right singular vector (state space), `u_min`
the left (residual space). Computed by Lanczos on `(JᵀJ)⁻¹` against a single
KLU factorization of `J`."""
function _min_singular_triplet(
    data::ACPowerFlowData,
    time_step::Int;
    x0::Union{Vector{Float64}, Nothing} = nothing,
    tol::Float64 = 1e-8,
    maxiter::Int = 300,
    krylovdim::Int = 40,
)
    residual = ACPowerFlowResidual(data, time_step)
    jac = ACPowerFlowJacobian(residual, time_step)
    x = isnothing(x0) ? calculate_x0(data, time_step) : copy(x0)
    residual(x, time_step)
    jac(time_step)

    cache = make_linear_solver_cache(PNM.KLUSolver(), jac.Jv)
    full_factor!(cache, jac.Jv)
    n = size(jac.Jv, 1)
    op = GramInverseOperator(cache, Vector{Float64}(undef, n))
    v0 = fill(1.0 / sqrt(n), n)
    vals, vecs, _ = KrylovKit.eigsolve(
        op, v0, 1, :LM; tol = tol, maxiter = maxiter, krylovdim = krylovdim,
        issymmetric = true,
    )
    σ_min = 1.0 / sqrt(abs(vals[1]))          # 1/σ_min² is the dominant eigenvalue
    v_min = vecs[1] ./ norm(vecs[1])          # right singular vector (state space)
    u = jac.Jv * v_min
    u_min = u ./ norm(u)                      # left singular vector (residual space)
    # F at the same state the mode was taken at. Snapshotted because `residual.Rv` is
    # scratch that any later residual evaluation overwrites.
    F = copy(residual.Rv)
    # ‖J⁻¹F‖: the full Newton step, from the same factorization. Used to compute
    # "how much of the step is the near-null direction"
    step = copy(F)
    newton_step_norm = try
        solve!(cache, step)
        norm(step)
    catch
        NaN   # singular J: the step is undefined, not zero
    end
    return (; σ_min, v_min, u_min, F, newton_step_norm, residual, jac)
end

"""
    find_jacobian_null_space(data, time_step; x0, k, tol, maxiter, krylovdim)
        -> (σ_min, state_entries, residual_entries)

Locate the near-null space of the AC power flow Jacobian `J` at `x0` (default:
the flat start). Computes the smallest singular value `σ_min` and its singular
vectors via Lanczos on `(JᵀJ)⁻¹` (dominant eigenvalue `1/σ_min²`), reusing a KLU
factorization of `J`. A tiny `σ_min` is the true rank-deficiency behind a huge
condition estimate `κ̂ ≈ σ_max/σ_min`; the returned vectors say *where* it lives.

Returns:
- `σ_min`: smallest singular value of `J`.
- `state_entries`: the `k` largest-magnitude entries of the right singular vector
  `v_min` (state directions `J` barely moves) as `(bus label, weight)` pairs.
- `residual_entries`: the `k` largest of the left singular vector
  `u_min = J·v_min / σ_min` (residual equations that are nearly unreachable), as
  `(formulation-aware label, weight)` pairs.
"""
function find_jacobian_null_space(
    data::ACPowerFlowData,
    time_step::Int;
    x0::Union{Vector{Float64}, Nothing} = nothing,
    k::Int = 10,
    tol::Float64 = 1e-8,
    maxiter::Int = 300,
    krylovdim::Int = 40,
)
    (; σ_min, v_min, u_min, residual) =
        _min_singular_triplet(data, time_step; x0, tol, maxiter, krylovdim)
    state_entries = _top_entries(v_min, k, ix -> _state_bus_label(residual, data, ix))
    residual_entries =
        _top_entries(
            u_min,
            k,
            ix -> _describe_residual_entry(residual, data, time_step, ix),
        )
    return σ_min, state_entries, residual_entries
end

# ---------------------------------------------------------------------------
# Bottleneck localizer: turn the critical mode into a cutset + reactive-reserve
# statement. Baseline-free — the weak boundary is a structural property of the
# operating region, not of whatever perturbation exposed it.
# ---------------------------------------------------------------------------

# original-bus-number → surviving (reduced) bus number, from the forward
# `bus_reduction_map: survivor → Set{absorbed}`. Buses that survive map to
# themselves. Trivial (no reduction) ⇒ empty dict; callers fall back to identity.
function _orig_to_survivor(data::ACPowerFlowData)
    reduction = get_network_reduction_data(data).bus_reduction_map
    o2s = Dict{Int, Int}()
    for (survivor, absorbed) in reduction
        o2s[survivor] = survivor
        for a in absorbed
            o2s[a] = survivor
        end
    end
    return o2s
end

# Reactive-power limits as `(min, max)`, or `nothing` if the device has none
# (e.g. a renewable/source without Q limits). Tolerant across device types.
function _gen_q_limits(g)
    try
        lims = PSY.get_reactive_power_limits(g)
        return lims === nothing ? nothing : lims
    catch
        return nothing
    end
end

"""
    _arc_flow_sensitivity(data, Δvm, Δθ, time_step) -> Vector{Float64}

Per-arc `|ΔS|`: the first-order change in branch complex power flow when the bus state
moves along the mode `(Δvm, Δθ)`. Exact analytic derivative, not a finite difference.

With `V = vm·e^{iθ}` the perturbation is `dV = (Δvm + i·vm·Δθ)·e^{iθ}`, and from the
package's own flow convention (`Sft = V_f · conj(Yft·V)`, cf. `solve_ac_power_flow.jl`):

    dSft = dV_f · conj(Yft·V) + V_f · conj(Yft·dV)

Both ends are computed and the larger magnitude is kept: a tear shows at both ends, and
they differ by the branch's own loss/charging terms, so taking the max avoids
under-ranking a branch measured at its lighter end.

Cost is four sparse mat-vecs against matrices that already exist. The operating point is
read from `data.bus_magnitude`/`bus_angles`, which the preceding residual evaluation in
[`_min_singular_triplet`](@ref) has already synced to the state the mode was taken at.
"""
function _arc_flow_sensitivity(
    data::ACPowerFlowData,
    Δvm::Vector{Float64},
    Δθ::Vector{Float64},
    time_step::Int,
)
    Yft = data.power_network_matrix.arc_admittance_from_to
    Ytf = data.power_network_matrix.arc_admittance_to_from
    arcs = PNM.get_arc_axis(Yft)
    bus_lookup = get_bus_lookup(data)
    fb_ix = [bus_lookup[first(a)] for a in arcs]
    tb_ix = [bus_lookup[last(a)] for a in arcs]

    vm = @view data.bus_magnitude[:, time_step]
    θ = @view data.bus_angles[:, time_step]
    rot = cis.(θ)
    V = vm .* rot
    dV = (Δvm .+ im .* vm .* Δθ) .* rot

    # PNM stores these as ComplexF32; widen so the sensitivity isn't quantized to
    # single precision (cf. the ~1e-4 flow discrepancy noted in `solve_ac_power_flow.jl`).
    Yft_d = ComplexF64.(Yft.data)
    Ytf_d = ComplexF64.(Ytf.data)
    dSft = dV[fb_ix] .* conj.(Yft_d * V) .+ V[fb_ix] .* conj.(Yft_d * dV)
    dStf = dV[tb_ix] .* conj.(Ytf_d * V) .+ V[tb_ix] .* conj.(Ytf_d * dV)
    return max.(abs.(dSft), abs.(dStf))
end

const _LCC_ROW_SYMBOLS = (:P_setpoint, :dc_balance, :rectifier_alpha, :inverter_alpha)

"""Structured (machine-parseable) classification of a polar residual entry as
`(kind, id, quantity)`: `(:bus, bus_number, :P|:Q)` for the bus block, or
`(:lcc, (from, to), row_symbol)` for the LCC tail. The companion to
[`_describe_residual_entry`](@ref), which returns the human-readable string."""
function _classify_residual_entry(
    ::ACPowerFlowResidual,
    data::ACPowerFlowData,
    ::Int,
    ix::Int,
)
    n_bus_eqs = 2 * size(data.bus_type, 1)
    if ix <= n_bus_eqs
        bus_ix = div(ix - 1, 2) + 1
        return (:bus, _diag_bus_number(data, bus_ix), isodd(ix) ? :P : :Q)
    end
    tail = ix - n_bus_eqs
    i = div(tail - 1, 4) + 1
    return (:lcc, data.lcc.arcs[i], _LCC_ROW_SYMBOLS[mod1(tail, 4)])
end

"""How many leading entries of a descending-sorted magnitude vector are
"significant" — i.e. within a factor `significance` of the peak. This replaces a
hard top-`k` with a data-driven cut: a lone dominant entry (the next is `1/200`
the peak) yields 1, a flat profile (many comparable entries) yields many. Always
returns at least 1 and never more than `k_max`."""
function _significant_count(sorted_desc::AbstractVector{<:Real};
    significance::Float64, k_max::Int)
    isempty(sorted_desc) && return 0
    peak = float(sorted_desc[1])
    peak <= 0 && return 1
    n = 0
    for v in sorted_desc
        v >= significance * peak || break
        n += 1
        n >= k_max && break
    end
    return max(n, 1)
end

# TODO: docstring out-of-date.
"""
    localize_bottleneck(data, sys, time_step; x0, k, q_margin_frac, verbose)
        -> (; σ_min, alignment, mode_step, mode_step_fraction, newton_step_norm,
             pocket, cutset, binding, exhausted_q, cuts)

`mode_step = |u_minᵀF|/σ_min`: the magnitude of the critical mode's contribution to the Newton step.

`mode_step_fraction = mode_step / ‖J⁻¹F‖ ∈ [0, 1]`: how much of the Newton step is that mode.
Values near 1 indicate that a fold is blowing up the step.

`alignment = |u_minᵀF|/‖F‖  ∈ [0, 1]`: how much of the residual is unreachable. Values near
1 indicate that the residual vector is aligned with the near-null space of `J`. In contrast,
if `mode_step_fraction` is near 1 but `alignment` is near 0, then the singular-ness of `J`
(very large `1/σ_min`) is to blame, not the residual's direction.

Turn the near-singular critical mode of the AC power flow Jacobian into a
*systemic bottleneck* statement. Fairly efficient: reuses a single KLU factorization.

Adaptive cutting: an entry is kept while it stays within a factor `significance` (default 0.1) 
of the list's peak, capped at `k`. A lone dominant bus/branch yields a one-element list;
a flat profile yields several. See [`_significant_count`](@ref).

Polar formulation only. Pinned quantities (e.g. `theta` at REF) mode component of exactly 0,
and the generator-power slots are dropped from the pocket entirely. Returns:
- `pocket`: `(bus_number, ΔVm, Δθ, driver)` the problematic state variables, read off
   from the right singular vector `v_min`. (We skip the `ΔP` and `ΔQ` state variables.)
   `|ΔVm|` and `|Δθ|` are cut independently and reported set is their union.
   `driver ∈ (:V, :θ, :both)` says which cut admitted the bus.
   Rows are ordered by `max(|ΔVm|/peak_V, |Δθ|/peak_θ)`.
- `cutset`: `(branch_name, from, to, |ΔS|, Δθ_mode, x)`: problematic branches, as evidenced
   by high first-order change in power flow in the direction of `v_min` . `Δθ_mode` and `x`
   are retained as context.
- `binding`: the problematic power-balance equations, read off from the left singular vector `u_min`,
   as machine-parseable `(kind, id, quantity, weight, F, contribution)`. Here
  `kind/id/quantity` is either `(:bus, bus_number, :P|:Q)` or `(:lcc, (from, to), row)`.
   Ranked and cut on `contribution = u_min[i]·F[i]`.  Large `weight` alone flags equations
  that are structurally hard even where no mismatch sits; large `F` alone flags mismatch the
  solver may clear easily. [TODO: when `F` is zero or near-zero, `u_min[i]` would be a more
  appropriate ranking metric.]
- `exhausted_q`: `(name, bus, Q, Qmin, Qmax, at_limit, via)` the problematic reactive power
  buses. Here, buses meet 2 criteria. First, they must be at or near one of their reactive
  power limits. Second, they must either have large `max(|ΔVm|,|Δθ|)` in `v_min` (via is `:pocket`),
  or large `|ΔQ_gen|` in `v_min` (via is `:q_mode`).
- `verbose` (default `true`) also `@info`-logs a readable summary.
"""
function localize_bottleneck(
    data::ACPowerFlowData,
    sys::PSY.System,
    time_step::Int = 1;
    x0::Union{Vector{Float64}, Nothing} = nothing,
    significance::Float64 = 0.1,
    k::Int = 15,
    binding_k::Int = 100,
    q_margin_frac::Float64 = 0.02,
    q_range_tol::Float64 = 1e-6,
    verbose::Bool = true,
)
    (; σ_min, v_min, u_min, F, newton_step_norm, residual) =
        _min_singular_triplet(data, time_step; x0)
    residual isa ACPowerFlowResidual || error(
        "localize_bottleneck supports the polar formulation only; got $(typeof(residual)).",
    )

    n_bus = size(data.bus_type, 1)
    bus_types = view(data.bus_type, :, time_step)
    # drop state variables that are generation; use 0.0 as Δθ at ref bus since angle is fixed.
    Δvm = [bus_types[i] == PSY.ACBusTypes.PQ ? v_min[2i - 1] : 0.0 for i in 1:n_bus]
    Δθ = [bus_types[i] == PSY.ACBusTypes.REF ? 0.0 : v_min[2i] for i in 1:n_bus]
    θ_mode = Δθ
    # Generator reactive participation, for the reactive power scan below.
    ΔQgen = [
        if bus_types[i] == PSY.ACBusTypes.PV
            v_min[2i - 1]
        elseif bus_types[i] == PSY.ACBusTypes.REF
            v_min[2i]
        else
            0.0
        end
        for i in 1:n_bus
    ]

    bus_numbers = collect(axes(data.power_network_matrix, 1))   # reduced ix → number
    num_to_ix = Dict(no => i for (i, no) in enumerate(bus_numbers))

    # Pocket: cut |ΔVm| and |Δθ| independently, then report the union.
    k_pocket = min(k, n_bus)
    peak_v, peak_θ = maximum(abs, Δvm), maximum(abs, Δθ)
    peak_q = maximum(abs, ΔQgen)
    ord_v = sortperm(Δvm; by = abs, rev = true)
    ord_θ = sortperm(Δθ; by = abs, rev = true)
    ord_q = sortperm(ΔQgen; by = abs, rev = true)
    mag_v, mag_θ = abs.(Δvm[ord_v]), abs.(Δθ[ord_θ])
    mag_q = abs.(ΔQgen[ord_q])
    # `_significant_count` floors at 1, which would admit an all-zero coordinate; an
    # identically-zero mode component has no participating buses at all.
    nv = peak_v > 0 ? _significant_count(mag_v; significance, k_max = k_pocket) : 0
    nθ = peak_θ > 0 ? _significant_count(mag_θ; significance, k_max = k_pocket) : 0
    set_v, set_θ = Set(ord_v[1:nv]), Set(ord_θ[1:nθ])
    # Rank-normalize so the two coordinates are comparable for display reasons.
    # No impact on cuts: decides which row prints first.
    _score(i) = max(
        peak_v > 0 ? abs(Δvm[i]) / peak_v : 0.0,
        peak_θ > 0 ? abs(Δθ[i]) / peak_θ : 0.0,
    )
    pocket_ix = sort!(collect(union(set_v, set_θ)); by = _score, rev = true)
    _driver(i) = i in set_v ? (i in set_θ ? :both : :V) : :θ
    pocket = [
        (bus_numbers[i], round(Δvm[i]; sigdigits = 4), round(Δθ[i]; sigdigits = 4),
            _driver(i))
        for i in pocket_ix
    ]

    # Cutset: branches ranked by |ΔS|, the first-order change in branch complex power
    # flow along the mode.
    dS_arc = _arc_flow_sensitivity(data, Δvm, Δθ, time_step)
    arc_ix = Dict(
        a => i for (i, a) in enumerate(
            PNM.get_arc_axis(data.power_network_matrix.arc_admittance_from_to))
    )

    o2s = _orig_to_survivor(data)
    surv_ix(no) = get(num_to_ix, get(o2s, no, no), 0)
    cut = NamedTuple{
        (:name, :from, :to, :dS, :dθ, :x),
        Tuple{String, Int, Int, Float64, Float64, Float64},
    }[]
    # TODO: better incorporate network reductions. Iterate over arcs, not system branches.
    for b in PSY.get_components(PSY.ACBranch, sys)
        occursin("ThreeWinding", String(nameof(typeof(b)))) && continue   # no single arc
        applicable(PSY.get_arc, b) || continue
        arc = PSY.get_arc(b)
        f0 = PSY.get_number(PSY.get_from(arc))
        t0 = PSY.get_number(PSY.get_to(arc))
        fi, ti = surv_ix(f0), surv_ix(t0)
        (fi == 0 || ti == 0 || fi == ti) && continue   # off-network or intra-cluster
        # Arcs are keyed by REDUCED bus numbers, in either orientation. Parallel
        # branches share one arc, so they share a |ΔS| — the aggregate is the right
        # quantity for a cutset, and each branch keeps its own name for actionability.
        key = (bus_numbers[fi], bus_numbers[ti])
        a = get(arc_ix, key, get(arc_ix, reverse(key), 0))
        a == 0 && continue
        push!(
            cut,
            (
                name = PSY.get_name(b), from = f0, to = t0,
                dS = dS_arc[a], dθ = abs(θ_mode[fi] - θ_mode[ti]),
                x = applicable(PSY.get_x, b) ? PSY.get_x(b) : NaN,
            ),
        )
    end
    sort!(cut; by = c -> c.dS, rev = true)
    nc = _significant_count([c.dS for c in cut]; significance, k_max = k)
    cutset = [
        (c.name, c.from, c.to, round(c.dS; sigdigits = 4),
            round(c.dθ; sigdigits = 4), round(c.x; sigdigits = 4))
        for c in cut[1:min(nc, length(cut))]
    ]

    # Binding mismatch: which power-balance equations are nearly unreachable. Ranked and
    # cut on contribution cᵢ = u_min[i]·F[i]. The Newton step is Δx = −Σᵢ(uᵢᵀF/σᵢ)vᵢ,
    # so what blows up is Σᵢcᵢ = u_minᵀF amplified by 1/σ_min.  Note that
    # cᵢ is signed and so terms can cancel.
    contribution = u_min .* F
    nb = _significant_count(
        sort(abs.(contribution); rev = true); significance, k_max = binding_k)
    binding = map(partialsortperm(contribution, 1:nb; by = abs, rev = true)) do ix
        kind, id, quantity = _classify_residual_entry(residual, data, time_step, ix)
        (kind, id, quantity, round(u_min[ix]; sigdigits = 4),
            round(F[ix]; sigdigits = 4),
            round(contribution[ix]; sigdigits = 4))
    end
    # The critical mode's contribution to the Newton step, ‖Δx_min‖ = |u_minᵀF|/σ_min,
    # and its share of the full step. 
    F_norm = norm(F)
    alignment = F_norm > 0 ? abs(LinearAlgebra.dot(u_min, F)) / F_norm : NaN
    mode_step = σ_min > 0 ? abs(LinearAlgebra.dot(u_min, F)) / σ_min : NaN
    mode_step_fraction = if isfinite(newton_step_norm) && newton_step_norm > 0
        mode_step / newton_step_norm
    else
        NaN
    end

    # Reactive reserve. Cut |ΔQ_gen| on its own and scan the union with the pocket. 
    nq = peak_q > 0 ? _significant_count(mag_q; significance, k_max = k_pocket) : 0
    reduction = get_network_reduction_data(data).bus_reduction_map
    # Expand a surviving-bus set to the originals reduced into it, so generators sitting
    # on absorbed buses are still caught.
    function _expand(ixs)
        s = Set(bus_numbers[i] for i in ixs)
        for b in collect(s)
            haskey(reduction, b) && union!(s, reduction[b])
        end
        return s
    end
    pocket_numbers = _expand(pocket_ix)
    q_mode_numbers = _expand(ord_q[1:nq])
    scan_numbers = union(pocket_numbers, q_mode_numbers)

    exhausted_q = Tuple{String, Int, Float64, Float64, Float64, Symbol, Symbol}[]
    for g in PSY.get_components(PSY.Generator, sys)
        PSY.get_available(g) || continue
        g_bus = PSY.get_number(PSY.get_bus(g))
        g_bus in scan_numbers || continue
        lims = _gen_q_limits(g)
        lims === nothing && continue
        # skip generators with no reactive range at all
        lims.max - lims.min > q_range_tol || continue
        q = PSY.get_reactive_power(g)
        hi = lims.max - q
        lo = q - lims.min
        margin = max(abs(lims.max), abs(lims.min), 1.0) * q_margin_frac
        (hi <= margin || lo <= margin) || continue
        in_p, in_q = g_bus in pocket_numbers, g_bus in q_mode_numbers
        push!(
            exhausted_q,
            (PSY.get_name(g), g_bus, q, lims.min, lims.max,
                hi <= margin ? :max : :min,
                in_p ? (in_q ? :both : :pocket) : :q_mode),
        )
    end

    # Carry the cut parameters and the backstop verdicts with the result, so a consumer
    # can tell if we hit `k`/`binding_k` while still on a plateau. The pocket flags OR
    # together: either coordinate's cut hitting the backstop undercounts the union.
    truncated = (;
        pocket = _backstop_truncated(mag_v, nv, k_pocket, significance) ||
                 _backstop_truncated(mag_θ, nθ, k_pocket, significance),
        cutset = _backstop_truncated([c.dS for c in cut], nc, k, significance),
        binding = _backstop_truncated(
            sort(abs.(contribution); rev = true), nb, binding_k, significance),
        # Not a reported list but a SCAN GATE: truncating it silently narrows which
        # generators are examined, so a missed pinned unit looks like "no Q problem".
        q_mode = _backstop_truncated(mag_q, nq, k_pocket, significance),
    )
    cuts = (;
        significance, k, binding_k, truncated,
        n_pocket_candidates = n_bus,
        n_cutset_candidates = length(cut),
        n_binding_candidates = length(u_min),
        n_q_mode_candidates = n_bus,
    )
    result = (; σ_min, alignment, mode_step, mode_step_fraction, newton_step_norm,
        pocket, cutset, binding, exhausted_q, cuts)
    verbose && _report_bottleneck(result)
    return result
end

"""
    _backstop_truncated(sorted_desc, n_kept, k_max, significance) -> Bool

Whether a [`_significant_count`](@ref) cut of `sorted_desc` that kept `n_kept` entries was 
made by the runaway *backstop* rather than by the data.

True iff all three hold: `n_kept == k_max`, the pool held more than `k_max` candidates
(so there was something left to cut), and the last kept magnitude is still within
`significance` of the peak (so the data-driven cut never fired).

Takes the full candidate vector rather than the kept slice so the "was there more?"
test is exact — it is called at the cut site, where that vector is still in hand.
"""
function _backstop_truncated(
    sorted_desc::AbstractVector{<:Real},
    n_kept::Int,
    k_max::Int,
    significance::Float64,
)
    n_kept == k_max || return false
    n_kept >= 1 || return false
    k_max < length(sorted_desc) || return false
    peak = float(sorted_desc[1])
    peak > 0 || return false
    return float(sorted_desc[n_kept]) >= significance * peak
end

"""
    bottleneck_tables(b) -> (; pocket, cutset, binding, binding_summary, exhausted_q,
                              truncated)

Reshape a [`localize_bottleneck`](@ref) result into `DataFrame`s, one per section.
"""
function bottleneck_tables(b)
    pocket = DataFrames.DataFrame(;
        bus = [p[1] for p in b.pocket],
        dVm = [p[2] for p in b.pocket],
        dtheta = [p[3] for p in b.pocket],
        driver = [p[4] for p in b.pocket],
    )
    cutset = DataFrames.DataFrame(;
        branch = [c[1] for c in b.cutset],
        from = [c[2] for c in b.cutset],
        to = [c[3] for c in b.cutset],
        dS = [c[4] for c in b.cutset],
        dtheta_mode = [c[5] for c in b.cutset],
        x = [c[6] for c in b.cutset],
    )
    # `id` is heterogeneous by design — an `Int` bus number for `:bus` rows, a
    # `(from, to)` arc tuple for `:lcc` rows — so the column lands as `Any`. Left
    # raw rather than stringified so the frame stays filterable, not just printable.
    binding = DataFrames.DataFrame(;
        kind = [e[1] for e in b.binding],
        id = Any[e[2] for e in b.binding],
        quantity = [e[3] for e in b.binding],
        weight = [e[4] for e in b.binding],
        F = [e[5] for e in b.binding],
        contribution = [e[6] for e in b.binding],
    )
    binding_summary = if DataFrames.nrow(binding) == 0
        DataFrames.DataFrame(;
            kind = Symbol[], quantity = Symbol[], n = Int[],
            c_max = Float64[], c_min = Float64[], F_max = Float64[],
        )
    else
        s = DataFrames.combine(
            DataFrames.groupby(binding, [:kind, :quantity]),
            DataFrames.nrow => :n,
            # Ranges of the RANKING quantity, so the summary and the row order agree.
            :contribution => (c -> _sf4(maximum(abs, c))) => :c_max,
            :contribution => (c -> _sf4(minimum(abs, c))) => :c_min,
            :F => (f -> _sf4(maximum(abs, f))) => :F_max,
        )
        DataFrames.sort!(s, :c_max; rev = true)
    end
    exhausted_q = DataFrames.DataFrame(;
        generator = [e[1] for e in b.exhausted_q],
        bus = [e[2] for e in b.exhausted_q],
        q = [e[3] for e in b.exhausted_q],
        q_min = [e[4] for e in b.exhausted_q],
        q_max = [e[5] for e in b.exhausted_q],
        at_limit = [e[6] for e in b.exhausted_q],
        via = [e[7] for e in b.exhausted_q],
    )
    # Passed through from the cut site rather than recomputed — see `localize_bottleneck`.
    truncated = if hasproperty(b, :cuts)
        b.cuts.truncated
    else
        (; pocket = false, cutset = false, binding = false, q_mode = false)
    end
    return (; pocket, cutset, binding, binding_summary, exhausted_q, truncated)
end

"""Print `df` under `title`, capped at `max_rows` with an elided-row note. The cap
is a readability guard, not a silent truncation — the omitted count is always
stated, so a plateau that got clipped announces itself."""
function _show_bottleneck_table(io::IO, title::AbstractString, df; max_rows::Int)
    n = DataFrames.nrow(df)
    println(io, "\n", title, " (", n, " row", n == 1 ? "" : "s", ")")
    n == 0 && return
    shown = min(n, max_rows)
    # `summary = false` drops DataFrames' own "N×M DataFrame" line: the count is
    # already in the title, and there it can state the elision too.
    show(io, DataFrames.first(df, shown); summary = false, eltypes = false,
        allrows = true, truncate = 40)
    println(io)
    shown < n && println(io, "  … ", n - shown, " further row(s) not shown.")
    return
end

"""
    _report_bottleneck(b; io, max_rows)

Print a readable summary of a [`localize_bottleneck`](@ref) result as tables (see
[`bottleneck_tables`](@ref)). Called automatically by `localize_bottleneck` when
`verbose = true`.

`max_rows` (default 20) caps each table independently; anything elided is reported
as a count. The binding section prints its `(kind, quantity)` roll-up first, so a
100-row plateau reads as a few summary lines before the detail rows.
"""
function _report_bottleneck(b; io::IO = stdout, max_rows::Int = 20)
    t = bottleneck_tables(b)
    println(io, "="^72)
    println(io, "Bottleneck localization: σ_min(J) = ", _sf4(b.σ_min))
    if hasproperty(b, :mode_step)
        f = b.mode_step_fraction
        verdict = if isnan(f)
            "step undefined (singular J or zero residual)"
        elseif f >= 0.9
            "the Newton step is essentially PURE critical mode"
        elseif f >= 0.5
            "the critical mode DOMINATES the Newton step"
        else
            "$(_sf4(100f))% of the step"
        end
        println(io, "  critical-mode step |u_minᵀF|/σ_min = ", _sf4(b.mode_step),
            "  of  ‖J⁻¹F‖ = ", _sf4(b.newton_step_norm), "   (", verdict, ")")
        println(io, "  alignment |u_minᵀF|/‖F‖ = ", _sf4(b.alignment),
            "  (unreachable fraction of the RESIDUAL)")
    elseif hasproperty(b, :alignment)
        println(io, "  alignment |u_minᵀF|/‖F‖ = ", _sf4(b.alignment))
    end
    println(io, "="^72)
    _show_bottleneck_table(io, "Collapse pocket (buses that collapse together)",
        t.pocket; max_rows)
    _show_bottleneck_table(io, "Weak cutset (branches the mode tears across)",
        t.cutset; max_rows)
    _show_bottleneck_table(io, "Binding mismatch by category",
        t.binding_summary; max_rows)
    _show_bottleneck_table(io, "Binding mismatch (nearly-infeasible equations)",
        t.binding; max_rows)
    if DataFrames.nrow(t.exhausted_q) == 0
        println(io,
            "\nNo pocket generators near a Q-limit " *
            "(bottleneck is transfer/structural, not local Q).")
    else
        _show_bottleneck_table(io, "Pocket generators out of reactive reserve",
            t.exhausted_q; max_rows)
    end
    println(io, "="^72)
    _warn_backstop_truncation(b, t)
    return
end

"""`@warn` for each bottleneck list that hit its runaway backstop while still on a
plateau (see [`_backstop_truncated`](@ref)) — the region is larger than reported.
Deliberately `@warn` rather than part of the printed tables: an unattended pipeline
consumes `bottleneck_tables` and never reads `stdout`, and an undercounted region is
a wrong *answer*, not a display artifact."""
function _warn_backstop_truncation(b, t = bottleneck_tables(b))
    hasproperty(b, :cuts) || return
    c = b.cuts
    for (field, knob, limit) in (
        (:pocket, "k", c.k),
        (:cutset, "k", c.k),
        (:binding, "binding_k", c.binding_k),
        (:q_mode, "k", c.k),
    )
        hasproperty(t.truncated, field) || continue   # pre-`q_mode` result
        getproperty(t.truncated, field) || continue
        n_cand = getproperty(c, Symbol(:n_, field, :_candidates))
        what = if field === :q_mode
            "reactive-reserve scan gate (`ΔQ_gen`) was narrowed by"
        else
            "`$field` list was cut by"
        end
        tail = if field === :q_mode
            "so generators outside it were never examined and a \"no Q problem\" " *
            "verdict may be an artifact of the gate"
        else
            "so the region continues past what is reported and is undercounted"
        end
        @warn "localize_bottleneck: the $what its `$knob = $limit` backstop, not by " *
              "the data — every entry kept is still within " *
              "`significance = $(c.significance)` of the peak, $tail. Re-run with a " *
              "larger `$knob` to see where it ends."
    end
    return
end

"""Run one iteration's diagnostics against the current `J`/residual. Does the
*single* per-iteration refactor of `cache` on `J.Jv` (NR/TR pass `linSolveCache`, LM
its own KLU `diag_cache`) and, on success, the *single* eigensolve shared by the
monitor line (`monitor`) and the fold bail-out (`bail`). Returns `true` iff the
caller should abort. A `SingularException` is itself a fold signature: under `bail`
it aborts, under monitor-only it reports `singular` and continues; any other
exception is rethrown."""
function run_solver_diagnostics!(
    state::SolverDiagnosticsState,
    label::AbstractString,
    residual::Union{ACPowerFlowResidual, ACRectangularCIResidual, ACMixedCPBResidual},
    J::Union{ACPowerFlowJacobian, ACRectangularCIJacobian, ACMixedCPBJacobian},
    time_step::Int,
    cache::PFLinearSolverCache,
    monitor::Bool,
    bail::Bool,
)::Bool
    # KLU throws `SingularException` on a singular J; AppleAccelerate does not (it
    # silently returns garbage but still factors), so only KLU reaches the catch.
    singular = false
    try
        numeric_refactor!(cache, J.Jv)
    catch e
        e isa LinearAlgebra.SingularException || rethrow()
        singular = true
    end

    data = J.data
    if singular
        if bail
            @warn "$label: the Jacobian is singular; this is a fold / " *
                  "voltage-collapse signature, aborting the search."
            return true
        end
        # Monitor-only: report the singularity rather than crashing, and leave
        # `state.prev_F` untouched so the next contraction ratio is meaningful.
        abs_max, ix = findmax(abs, residual.Rv)
        @info "$label: ‖F‖_∞ = $(_sf4(abs_max)) at " *
              "$(_describe_residual_entry(residual, data, time_step, ix)), " *
              "κ̂(J) = singular, λ_min(S) = singular"
        return false
    end

    n_state = size(J.Jv, 1)
    # Trailing block is the FULL non-bus tail (LCC + VSC + area interchange), not just
    # LCC: n_state on a VSC/area-interchange system is larger than 2*nbuses + 4*n_lcc.
    n_bus = n_state - state_tail_length(data, get_dc_network(data))
    op = SchurInverseOperator(cache, n_bus, state.buffer)
    λ_min, converged = _schur_min_eigenvalue(op)

    if monitor
        abs_max, ix = findmax(abs, residual.Rv)
        κ = _diag_condest(cache)
        λ_str = if converged
            "$(_fmt_eig(λ_min)) (|λ_min| = $(_sf4(abs(λ_min))))"
        else
            "not-converged"
        end
        parts = [
            "‖F‖_∞ = $(_sf4(abs_max)) at " *
            "$(_describe_residual_entry(residual, data, time_step, ix))",
            "κ̂(J) = $(isnan(κ) ? "n/a (KLU-only)" : string(_sf4(κ)))",
            "λ_min(S) = $λ_str",
        ]
        if !isnan(state.prev_F) && state.prev_F > 0
            push!(parts, "contraction = $(_sf4(abs_max / state.prev_F))")
        end
        @info "$label: " * join(parts, ", ")
        state.prev_F = abs_max
    end

    if bail
        return _decide_eig_sign_switch!(state, label, λ_min, converged)
    end
    return false
end

"""
    _report_area_interchange_failure(data, time_step)

Terminal-failure diagnostic for embedded area net-interchange control, called when the
greedy relax loop exhausts the enrolled set without converging. The WORKING set is empty
at that point, so this reports against the PRISTINE tie/area set at the last attempted
iterate's bus state, naming the area with the largest-magnitude interchange-row residual.
No-op if area interchange control was never enrolled.
"""
function _report_area_interchange_failure(data::ACPowerFlowData, time_step::Int)
    aid = data.area_interchange
    isempty(aid.pristine_areas) && return
    gaps = [
        _area_net_interchange(
            aid.pristine_ties, aid.pristine_dc_ties, area.tail_ix, data, time_step,
        ) - area.pdes
        for area in aid.pristine_areas
    ]
    abs_max, ix = findmax(abs, gaps)
    area = aid.pristine_areas[ix]
    @warn "Area interchange: Newton did not converge after the greedy relax loop " *
          "de-enrolled every controlled area (network non-convergence, not a relaxed " *
          "schedule); the largest interchange-row residual at the last attempted " *
          "iterate is area \"$(area.name)\" with |r| = $(_sf4(abs_max)) " *
          "(target PDES = $(area.pdes))."
    return
end
