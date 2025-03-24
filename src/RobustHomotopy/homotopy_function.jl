struct HomotopyFunction
    data::ACPowerFlowData
    last_x::Vector{Float64}
    pfResidual::ACPowerFlowResidual
    t_k_ref::Base.RefValue{Float64}
end

function HomotopyFunction(data::ACPowerFlowData)
    time_step = 1
    pfResidual = ACPowerFlowResidual(data, time_step)
    t_k = 0.0
    return HomotopyFunction(data, similar(pfResidual.Rv), pfResidual, Ref(t_k))
end

function (homFunc::HomotopyFunction)(x::Vector{Float64})
    copyto!(homFunc.last_x, x)
    time_step = 1
    homFunc.pfResidual(x, time_step)
    t_k = homFunc.t_k_ref[]
    oddStateVars = @view x[1:2:end]
    PQ_mask = get_bus_type(homFunc.data)[:, time_step] .== (PSY.ACBusTypes.PQ, )
    return 0.5*(1 - t_k)*sum(abs2, oddStateVars[PQ_mask] .- 1.0)
                                + 0.5*t_k*sum(abs2, pfResidual.Rv)
end