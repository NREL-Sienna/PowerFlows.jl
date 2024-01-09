struct ContinuationPowerFlow{T <: PolarPowerFlowJacobian, U <: PolarPowerFlow}
    jacobian::T
    pf::U
    extended_J::SparseArrays.SparseMatrixCSC{Float64, Int32}
    data::ACPowerFlowData
    load_increase_delta::Float64
    change_factor::Float64
    lambda::Base.RefValue{Float64}
    tolerance::Float64
    maxiterations::Float64
    generator_participation_factors::Vector{Float64}
    
end

function ContinuationPowerFlow(
    sys::PSY.System;
    tolerance::Float64=1e-3,
    maxiterations::Int=25,
    change_factor::Float64,
    generator_participation_factors::Vector{Float64} = Float64[],
    check_reactive_power_limits = true,
    check_connectivity = true)

    data = PowerFlowData(
        ACPowerFlow(; check_reactive_power_limits = check_reactive_power_limits),
        system;
        check_connectivity = check_connectivity,
    )
    data = ACPowerFlowData(sys)
    pf = PolarPowerFlow(data)
    J = PowerFlows.PolarPowerFlowJacobian(data, pf.x0)

    return ContinuationPowerFlow(
        J,
        pf,
        data,

    )

end


function run!(cpf::ContinuationPowerFlow, bus::PSY.Bus)


end
