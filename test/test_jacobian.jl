@testset "Jacobian Verification" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg2000_sys")
    pf = PF.ACPowerFlow{NewtonRaphsonACPowerFlow}()
    data = PF.PowerFlowData(pf, sys; correct_bustypes = true)
    time_step = 1
    residual = PF.ACPowerFlowResidual(data, time_step)
    J = PF.ACPowerFlowJacobian(data, time_step)
    n = 2 * length(collect(get_components(ACBus, sys)))
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
