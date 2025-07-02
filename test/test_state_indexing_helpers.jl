@testset "partition state" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    pf = ACPowerFlow()
    data = PowerFlowData(pf, sys)
    time_step = 1
    x0 = collect(1.0:1.0:10.0) # just need 10 distinct Float64's
    bus_types = @view PF.get_bus_type(data)[:, time_step]
    (Ps_1, Qs_1, Vms_1, Vas_1) = (Vector{Float64}() for _ in 1:4)
    i = 1
    for bt in bus_types
        if bt == PSY.ACBusTypes.REF
            push!(Ps_1, x0[i])
            i += 1
            push!(Qs_1, x0[i])
            i += 1
            push!(Vms_1, NaN)
            push!(Vas_1, NaN)
        elseif bt == PSY.ACBusTypes.PV
            push!(Qs_1, x0[i])
            i += 1
            push!(Vas_1, x0[i])
            i += 1
            push!(Ps_1, NaN)
            push!(Vms_1, NaN)
        elseif bt == PSY.ACBusTypes.PQ
            push!(Vms_1, x0[i])
            i += 1
            push!(Vas_1, x0[i])
            i += 1
            push!(Ps_1, NaN)
            push!(Qs_1, NaN)
        end
    end
    tp = PF.partition_state(x0, bus_types)
    @test isequal(tp[:P], Ps_1)
    @test isequal(tp[:Q], Qs_1)
    @test isequal(tp[:Vm], Vms_1)
    @test isequal(tp[:Va], Vas_1)
end
