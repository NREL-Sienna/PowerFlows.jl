"""
Tests for TxStepping homotopy continuation power flow.

1. Jacobian numerical verification via finite differences.
2. Integration tests: convergence, voltage/power accuracy vs AC NR.
"""

# ===== Jacobian Verification =====

function verify_tx_stepping_jacobian(sys::PSY.System; lambda::Float64 = 0.0)
    pf = TxSteppingPowerFlow()
    data = PF.PowerFlowData(pf, sys)
    time_step = 1

    n = size(data.bus_type, 1)
    bus_types = @view data.bus_type[:, time_step]
    n_pv = count(==(PSY.ACBusTypes.PV), bus_types)
    dim = 2 * n + n_pv

    # Initialize state at a non-trivial point (flat start with Vset magnitudes)
    V0 = data.bus_magnitude[:, time_step] .* exp.((im,) .* data.bus_angles[:, time_step])
    q_g0 = zeros(Float64, n)
    state = PF.TxSteppingSolverState(V0, q_g0)

    PNM.update_y_lambda!(data.power_network_matrix, lambda, 1e4)

    ref_mag = PF._ref_mag_for_bus(data, time_step)
    residual = PF.TxSteppingResidual(data, state, zeros(Float64, dim), ref_mag)
    jacobian = PF._build_tx_stepping_jacobian(data, state, ref_mag, time_step)

    # Evaluate residual and Jacobian at the initial state
    PF.compute_residual!(residual, time_step, lambda)
    R0 = copy(residual.Rv)
    PF._update_jacobian!(
        jacobian.Jv, jacobian.rhs, state, data, ref_mag, time_step, lambda,
    )
    J0 = copy(jacobian.Jv)

    # Helper to perturb state and evaluate residual
    function eval_residual_at(state_base, delta, dim, n)
        V_pert = copy(state_base.V)
        q_g_pert = copy(state_base.q_g)
        # First 2n entries map to V (interleaved real/imag)
        dV = reinterpret(ComplexF64, delta[1:(2 * n)])
        V_pert .+= dV
        # Remaining entries map to q_g at PV buses
        pv_ix = 1
        for (ix, bt) in enumerate(bus_types)
            bt != PSY.ACBusTypes.PV && continue
            q_g_pert[ix] += delta[2 * n + pv_ix]
            pv_ix += 1
        end
        tmp_state = PF.TxSteppingSolverState(V_pert, q_g_pert)
        tmp_residual =
            PF.TxSteppingResidual(data, tmp_state, zeros(Float64, dim), ref_mag)
        PF.compute_residual!(tmp_residual, time_step, lambda)
        return copy(tmp_residual.Rv)
    end

    Δx_start, Δx_stop = 1e-5, 1e-8
    for j in 1:dim
        u = zeros(dim)
        u[j] = 1.0
        Δx_mag = Δx_start
        close_enough = false
        while !close_enough && Δx_mag >= Δx_stop
            R_plus = eval_residual_at(state, Δx_mag .* u, dim, n)
            ΔR = R_plus .- R0
            floating_point_issues = all(isapprox.(ΔR, 0.0; atol = eps(Float32)))
            if floating_point_issues
                break
            end
            ∂R_∂u_numerical = ΔR ./ Δx_mag
            ∂R_∂u_symbolic = J0 * u
            if !isapprox(∂R_∂u_numerical, ∂R_∂u_symbolic; rtol = 1e-4, atol = 1e-6)
                Δx_mag /= 10.0
            else
                close_enough = true
            end
        end
        @test close_enough
    end
end

@testset "TxStepping Jacobian verification" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    verify_tx_stepping_jacobian(sys)
end

@testset "TxStepping Jacobian verification (lambda=0.5)" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    verify_tx_stepping_jacobian(sys; lambda = 0.5)
end

# ===== Integration Tests =====

@testset "TxStepping convergence on c_sys5" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    pf_tx = TxSteppingPowerFlow()
    result_tx = solve_power_flow(pf_tx, sys)
    @test !ismissing(result_tx)
    @test haskey(result_tx, "bus_results")
end

@testset "TxStepping matches AC NR on c_sys5" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # Solve with AC NR (reference)
    pf_ac = ACPowerFlow{NewtonRaphsonACPowerFlow}()
    result_ac = solve_power_flow(pf_ac, sys)

    # Solve with TxStepping
    pf_tx = TxSteppingPowerFlow()
    result_tx = solve_power_flow(pf_tx, sys)
    @test !ismissing(result_tx)

    bus_ac = result_ac["bus_results"]
    bus_tx = result_tx["bus_results"]

    # Voltage magnitudes and angles should match tightly
    @test isapprox(bus_tx.Vm, bus_ac.Vm; atol = DIFF_INF_TOLERANCE)
    @test isapprox(bus_tx.θ, bus_ac.θ; atol = DIFF_INF_TOLERANCE)
    # Power comparisons: Ybus (used by AC NR) stores entries in ComplexF32 while
    # YbusSplit now uses ComplexF64, so TxStepping is actually more precise.
    # The ~1e-3 diff in REF bus power reflects Float32 truncation in the AC reference.
    @test isapprox(bus_tx.P_gen, bus_ac.P_gen; atol = 2e-3)
    @test isapprox(bus_tx.Q_gen, bus_ac.Q_gen; atol = 2e-3)
end

@testset "TxStepping convergence on c_sys14" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")

    pf_tx = TxSteppingPowerFlow()
    result_tx = solve_power_flow(pf_tx, sys)
    @test !ismissing(result_tx)
    @test haskey(result_tx, "bus_results")
end

@testset "TxStepping matches AC NR on c_sys14" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")

    pf_ac = ACPowerFlow{NewtonRaphsonACPowerFlow}()
    result_ac = solve_power_flow(pf_ac, sys)

    pf_tx = TxSteppingPowerFlow()
    result_tx = solve_power_flow(pf_tx, sys)
    @test !ismissing(result_tx)

    bus_ac = result_ac["bus_results"]
    bus_tx = result_tx["bus_results"]

    @test isapprox(bus_tx.Vm, bus_ac.Vm; atol = DIFF_INF_TOLERANCE)
    @test isapprox(bus_tx.θ, bus_ac.θ; atol = DIFF_INF_TOLERANCE)
    @test isapprox(bus_tx.P_gen, bus_ac.P_gen; atol = 2e-3)
    @test isapprox(bus_tx.Q_gen, bus_ac.Q_gen; atol = 2e-3)
end

@testset "TxStepping solve_and_store_power_flow!" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # Get reference values from AC NR
    pf_ac = ACPowerFlow{NewtonRaphsonACPowerFlow}()
    sys_ac = deepcopy(sys)
    @test solve_and_store_power_flow!(pf_ac, sys_ac)

    # Solve with TxStepping and store
    pf_tx = TxSteppingPowerFlow()
    sys_tx = deepcopy(sys)
    @test solve_and_store_power_flow!(pf_tx, sys_tx)

    # Compare bus voltages stored in the system
    for bus in PSY.get_components(PSY.ACBus, sys_tx)
        name = PSY.get_name(bus)
        bus_ac = PSY.get_component(PSY.ACBus, sys_ac, name)
        @test isapprox(
            PSY.get_magnitude(bus),
            PSY.get_magnitude(bus_ac);
            atol = DIFF_INF_TOLERANCE,
        )
        @test isapprox(
            PSY.get_angle(bus),
            PSY.get_angle(bus_ac);
            atol = DIFF_INF_TOLERANCE,
        )
    end
end
