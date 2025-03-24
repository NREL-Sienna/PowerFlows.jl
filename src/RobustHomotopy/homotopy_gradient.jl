struct HomotopyGradient # I would call this HomotopyFunctionAndGradient but that's quite long.
    # the gradient requires knowing the residual, thus why I combine them.
    data::ACPowerFlowData
    pfResidual::ACPowerFlowResidual
    J::ACPowerFlowJacobian
    grad::Vector{Float64}
    t_k_ref::Base.RefValue{Float64}
end

function HomotopyGradient(data::ACPowerFlowData)
    time_step = 1
    pfResidual = ACPowerFlowResidual(data, time_step)
    J = ACPowerFlowJacobian(data, time_step)
    t_k = 0.0
    nbuses = size(get_bus_type(data), 1)
    return HomotopyGradient(data, pfResidual, J, zeros(2*nbuses), Ref(t_k))
end

function (homGrad::HomotopyGradient)(g::Vector{Float64},
    x::Vector{Float64};
    internal = false)
    println("gradient at $x")
    t_k = homGrad.t_k_ref[]
    time_step = 1
    homGrad.pfResidual(x, time_step)
    homGrad.J(time_step)
    g .= t_k * homGrad.J.Jv' * homGrad.pfResidual.Rv
    for (bus_ix, bt) in enumerate(get_bus_type(homGrad.data)[:, time_step])
        if bt == PSY.ACBusTypes.PQ
            g[2*bus_ix - 1] += (1.0 - t_k) * (x[2*bus_ix - 1] - 1.0)
        end
    end
    if !internal
        copyto!(homGrad.grad, g)
    end
    return nothing
end

function (homGrad::HomotopyGradient)(x::Vector{Float64})
    homGrad(homGrad.grad, x; internal = true)
end

