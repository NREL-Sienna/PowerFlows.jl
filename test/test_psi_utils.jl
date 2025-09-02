@testset "contributes active/reactive power" begin
    device_types = Type[PSY.StaticInjection]
    while !isempty(device_types)
        T = pop!(device_types)
        if !isabstracttype(T)
            instance = T(nothing)
            if PF.contributes_active_power(instance)
                @test hasmethod(PSY.get_active_power, Tuple{T})
                @test hasmethod(PF.active_power_contribution_type, Tuple{T})
            end
            # TODO awkward carve-out here, for FACTSControlDevice.
            # for those, reactive_power_required is the equivalent of get_reactive_power,
            # but it depends on control mode, and PSY isn't updated for those modes yet.
            if PF.contributes_reactive_power(instance) && T != PSY.FACTSControlDevice
                @test hasmethod(PSY.get_reactive_power, Tuple{T})
                @test hasmethod(PF.reactive_power_contribution_type, Tuple{T})
            end
        end
        append!(device_types, InteractiveUtils.subtypes(T))
    end
end
