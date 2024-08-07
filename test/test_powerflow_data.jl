@testset "PowerFlowData" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    PowerFlowData(ACPowerFlow(), sys)
    PowerFlowData(DCPowerFlow(), sys)
    PowerFlowData(PTDFDCPowerFlow(), sys)
    PowerFlowData(vPTDFDCPowerFlow(), sys)
end

@testset "PowerFlowData multiperiod" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    timesteps = 24
    # TODO: "multiperiod AC still to implement"
    PowerFlowData(DCPowerFlow(), sys; timesteps = timesteps)
    PowerFlowData(PTDFDCPowerFlow(), sys; timesteps = timesteps)
    PowerFlowData(vPTDFDCPowerFlow(), sys; timesteps = timesteps)
end

@testset "System <-> PowerFlowData round trip" begin
    # TODO currently only tested with ACPowerFlow
    # TODO test that update_system! errors if the PowerFlowData doesn't correspond to the system

    sys_original = build_system(PSISystems, "RTS_GMLC_DA_sys")
    data_original = PowerFlowData(ACPowerFlow(), sys_original)

    sys_modified = deepcopy(sys_original)
    modify_rts_system!(sys_modified)
    data_modified = PowerFlowData(ACPowerFlow(), sys_original)
    modify_rts_powerflow!(data_modified)

    # update_system! with unmodified PowerFlowData should result in system that yields unmodified PowerFlowData
    # (NOTE does NOT necessarily yield original system due to power redistribution)
    sys_null_updated = deepcopy(sys_original)
    PF.update_system!(sys_null_updated, data_original)
    data_null_updated = PowerFlowData(ACPowerFlow(), sys_null_updated)
    # TODO fix bug in `_reactive_power_redistribution_pv`, see https://github.com/NREL-Sienna/PowerFlows.jl/issues/44
    @test IS.compare_values(powerflow_match_fn, data_null_updated, data_original;
        exclude = Set([:bus_reactivepower_injection]))

    # Modified versions should not be the same as unmodified versions
    @test !@test_logs((:error, r"values do not match"),
        match_mode = :any, min_level = Logging.Error,
        IS.compare_values(powerflow_match_fn, data_original, data_modified))
    @test !@test_logs((:error, r"values do not match"),
        match_mode = :any, min_level = Logging.Error,
        IS.compare_values(powerflow_match_fn, sys_original, sys_modified))

    # Constructing PowerFlowData from modified system should result in data_modified
    @test IS.compare_values(
        powerflow_match_fn,
        PowerFlowData(ACPowerFlow(), sys_modified),
        data_modified,
    )

    # The big one: update_system! with modified PowerFlowData should result in sys_modified
    sys_modify_updated = deepcopy(sys_original)
    PF.update_system!(sys_modify_updated, data_modified)
    # TODO fix bug in `_reactive_power_redistribution_pv`, see https://github.com/NREL-Sienna/PowerFlows.jl/issues/44
    @test IS.compare_values(powerflow_match_fn, sys_modify_updated, sys_modified;
        exclude = Set([:reactive_power]))
end
