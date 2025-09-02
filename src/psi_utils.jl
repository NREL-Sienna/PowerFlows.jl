"""
Check if a device has attribute 'active_power' for active power consumption or generation.
"""
function contributes_active_power(::T) where {T <: PSY.Device}
    throw(
        IS.NotImplementedError(
            "contributes_active_power not implemented for this device type $(T)",
        ),
    )
end

IS.@scoped_enum(PowerContributionType, INJECTION = 1, WITHDRAWAL = 2)
function active_power_contribution_type(::T) where {T <: PSY.Device}
    throw(
        IS.NotImplementedError(
            "active_power_contribution_type not implemented for this device type $(T)",
        ),
    )
end

# most static components inject power
contributes_active_power(::PSY.StaticInjection) = true
active_power_contribution_type(::PSY.StaticInjection) = PowerContributionType.INJECTION
# carve out the exceptions: we handle constant impedance/current loads separately
contributes_active_power(::PSY.SynchronousCondenser) = false
contributes_active_power(::PSY.FixedAdmittance) = false
contributes_active_power(::PSY.SwitchedAdmittance) = false
contributes_active_power(::PSY.FACTSControlDevice) = false
# loads withdraw power.
active_power_contribution_type(::PSY.ElectricLoad) = PowerContributionType.WITHDRAWAL

function contributes_reactive_power(::T) where {T <: PSY.Device}
    throw(
        IS.NotImplementedError(
            "contributes_reactive_power not implemented for this device type $(T)",
        ),
    )
end

function reactive_power_contribution_type(::T) where {T <: PSY.Device}
    throw(
        IS.NotImplementedError(
            "reactive_power_contribution_sign not implemented for this device type $(T)",
        ),
    )
end

# most static components contribute reactive power
contributes_reactive_power(::PSY.StaticInjection) = true
reactive_power_contribution_type(::PSY.StaticInjection) = PowerContributionType.INJECTION
# handle constant impedance loads separately
contributes_reactive_power(::PSY.FixedAdmittance) = false
contributes_reactive_power(::PSY.SwitchedAdmittance) = false
# interconnecting converters do not support reactive power
contributes_reactive_power(::PSY.InterconnectingConverter) = false
# loads withdraw reactive power
reactive_power_contribution_type(::PSY.ElectricLoad) = PowerContributionType.WITHDRAWAL
