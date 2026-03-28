# PSS/e DCLF-style loss injection (POM Section 8.4)
# Pre-compute branch losses from the AC voltage profile and inject at sending-end buses.

# Helpers to extract tap ratio and phase shift from direct branches,
# following the same dispatch pattern as PNM._get_tap in Ybus.jl.
_branch_tap(::PSY.ACBranch) = 1.0
_branch_tap(br::PSY.TwoWindingTransformer) = PSY.get_tap(br)
_branch_tap(::PSY.Transformer2W) = 1.0

_branch_shift(::PSY.ACBranch) = 0.0
_branch_shift(br::PSY.TwoWindingTransformer) = PSY.get_α(br)

"""
    _get_arc_branch_params(data::ABAPowerFlowData) -> (rs, xs, taps, shifts)

Return vectors of series resistance, reactance, tap ratio, and phase shift
for each arc. Follows the same dispatch over NetworkReductionData branch maps
as [`_get_arc_resistances`](@ref).
"""
function _get_arc_branch_params(data::ABAPowerFlowData)
    nrd = get_network_reduction_data(data)
    arc_ax = get_arc_axis(data)
    n_arcs = length(arc_ax)
    rs = zeros(n_arcs)
    xs = zeros(n_arcs)
    taps = ones(n_arcs)
    shifts = zeros(n_arcs)
    for (ix, arc) in enumerate(arc_ax)
        branch = get(PNM.get_direct_branch_map(nrd), arc, nothing)
        if branch !== nothing
            rs[ix] = PSY.get_r(branch)
            xs[ix] = PSY.get_x(branch)
            taps[ix] = _branch_tap(branch)
            shifts[ix] = _branch_shift(branch)
            continue
        end
        parallel = get(PNM.get_parallel_branch_map(nrd), arc, nothing)
        if parallel !== nothing
            eq = PNM.get_equivalent_physical_branch_parameters(parallel)
            rs[ix] = PNM.get_equivalent_r(eq)
            xs[ix] = PNM.get_equivalent_x(eq)
            taps[ix] = PNM.get_equivalent_tap(eq)
            shifts[ix] = PNM.get_equivalent_shift(eq)
            continue
        end
        series = get(PNM.get_series_branch_map(nrd), arc, nothing)
        if series !== nothing
            eq = PNM.get_equivalent_physical_branch_parameters(series)
            rs[ix] = PNM.get_equivalent_r(eq)
            xs[ix] = PNM.get_equivalent_x(eq)
            taps[ix] = PNM.get_equivalent_tap(eq)
            shifts[ix] = PNM.get_equivalent_shift(eq)
            continue
        end
        winding = get(PNM.get_transformer3W_map(nrd), arc, nothing)
        if winding !== nothing
            rs[ix] = PNM.get_equivalent_r(winding)
            xs[ix] = PNM.get_equivalent_x(winding)
            taps[ix] = PNM.get_equivalent_tap(winding)
            shifts[ix] = PNM.get_equivalent_α(winding)
            continue
        end
        added = get(PNM.get_added_branch_map(nrd), arc, nothing)
        if added !== nothing
            z = 1 / added
            rs[ix] = real(z)
            xs[ix] = imag(z)
            continue
        end
        error("Arc $arc not found in any branch map.")
    end
    return rs, xs, taps, shifts
end

# No-op for non-ABA data types (AC, PTDF, vPTDF).
_populate_loss_injections!(::PowerFlowData, ::PSY.System) = nothing

"""
    _populate_loss_injections!(data::ABAPowerFlowData, sys::PSY.System)

Compute PSS/e DCLF-style loss injections from the AC voltage profile stored
in `sys` and write them into `data.initial_loss_injections`.

For each in-service branch k with from-bus i and to-bus j:

    g_k = r_k / (r_k² + x_k²)
    P_loss_k = g_k · (Vi²/tap² + Vj² − 2·Vi·Vj/tap · cos(θi − θj − shift))

Losses are withdrawn at the sending-end bus (the higher-angle side).
This is a single-pass, non-iterative computation matching PSS/e DCLF (POM §8.4).

When the system is flat-start (V=1, θ=0 everywhere), all losses are zero and
the method degenerates to standard lossless DCLF.
"""
function _populate_loss_injections!(data::ABAPowerFlowData, sys::PSY.System)
    loss_inj = data.initial_loss_injections
    loss_inj === nothing && return
    if data.loss_factors !== nothing
        @warn "Both loss_approximation_as_injection and calculate_loss_factors are enabled. " *
              "The pre-injected AC-derived losses use a different approximation than " *
              "the PTDF-based loss factors; results may be inconsistent."
    end
    loss_inj .= 0.0

    arc_ax = get_arc_axis(data)
    bus_lookup = get_bus_lookup(data)

    # Build bus_number → (V, θ) map from the system's current AC state.
    bus_voltage = Dict{Int, Tuple{Float64, Float64}}()
    for bus in PSY.get_components(PSY.ACBus, sys)
        bus_voltage[PSY.get_number(bus)] =
            (PSY.get_magnitude(bus), PSY.get_angle(bus))
    end

    rs, xs, taps, shifts = _get_arc_branch_params(data)

    for (ix, arc) in enumerate(arc_ax)
        from_bus_no, to_bus_no = arc
        r = rs[ix]
        x = xs[ix]
        tap = taps[ix]
        shift = shifts[ix]

        # Skip zero-impedance arcs.
        z2 = r^2 + x^2
        z2 == 0.0 && continue

        g = r / z2

        # Read AC voltages from the system (not from PowerFlowData, which forces V=1 for DC).
        Vi, θi = bus_voltage[from_bus_no]
        Vj, θj = bus_voltage[to_bus_no]

        δ = θi - θj - shift
        P_loss = g * (Vi^2 / tap^2 + Vj^2 - 2 * Vi * Vj / tap * cos(δ))

        # Sending end = higher-angle bus; default to from-bus on tie (PSS/e convention).
        from_ix = bus_lookup[from_bus_no]
        to_ix = bus_lookup[to_bus_no]
        sending_ix = (θi >= θj) ? from_ix : to_ix

        # Withdraw losses at sending end (negative injection).
        @views loss_inj[sending_ix, :] .-= P_loss
    end
    return
end
