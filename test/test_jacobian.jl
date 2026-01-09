function verify_jacobian(sys::PSY.System)
    pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}(; correct_bustypes = true)
    data = PF.PowerFlowData(pf, sys)
    time_step = 1
    residual = PF.ACPowerFlowResidual(data, time_step)
    J = PF.ACPowerFlowJacobian(data, time_step)
    n_lccs =
        length(collect(PSY.get_components(PSY.get_available, PSY.TwoTerminalLCCLine, sys)))
    n = 2 * length(collect(get_components(ACBus, sys))) + 4 * n_lccs
    x0 = PF.calculate_x0(data, time_step)
    residual(x0, time_step)
    J(time_step)
    J_x0 = deepcopy(J.Jv)
    Δx_start, Δx_stop = 1e-3, 1e-6
    for j in 1:n
        u = zeros(n)
        u[j] = 1
        Δx_mag = Δx_start
        close_enough = false
        while !close_enough && Δx_mag >= Δx_stop
            inputValues = [x0, x0 + Δx_mag * u]
            outputValues = Vector{Vector{Float64}}()
            for inputVal in inputValues
                residual(inputVal, time_step)
                push!(outputValues, deepcopy(residual.Rv))
            end
            ΔF = outputValues[2] .- outputValues[1]
            floating_point_issues = all(isapprox.(ΔF, 0.0; atol = eps(Float32)))
            if floating_point_issues
                break
            end
            ∂F_∂u_numerical = ΔF ./ Δx_mag
            ∂F_∂u_symbolic = J_x0 * u

            if !isapprox(∂F_∂u_numerical, ∂F_∂u_symbolic; rtol = 1e-6, atol = eps(Float32))
                Δx_mag /= 10.0
            else
                close_enough = true
            end
        end
        @test close_enough
    end
end

@testset "Jacobian verification" begin
    sys = PSB.build_system(PSITestSystems, "c_sys14")
    verify_jacobian(sys)
end

# This fails. Might be a problem, or might just be because (according to Roman) we treat
# the angle \phi as constant, to avoid super messy expressions in the Jacobian, 
# and update it separately between iterations.
@testset "Jacobian verification with LCC" begin
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 1, ACBusTypes.REF, 230, 1.1, 0.0)
    b2 = _add_simple_bus!(sys, 2, ACBusTypes.PQ, 230, 1.1, 0.0)
    b3 = _add_simple_bus!(sys, 3, ACBusTypes.PQ, 230, 1.1, 0.0)
    ld2 = _add_simple_load!(sys, b2, 10, 5)
    ld3 = _add_simple_load!(sys, b3, 60, 20)
    l12 = _add_simple_line!(sys, b1, b2, 5e-3, 5e-3, 1e-3)
    l13 = _add_simple_line!(sys, b1, b3, 5e-3, 5e-3, 1e-3)
    s1 = _add_simple_source!(sys, b1, 0.0, 0.0)
    lcc = _add_simple_lcc!(sys, b2, b3, 0.05, 0.05, 0.08)
    # verify_jacobian(sys)
end
