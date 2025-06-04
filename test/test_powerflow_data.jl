@testset "PowerFlowData" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    @test PowerFlowData(ACPowerFlow{LUACPowerFlow}(), sys; fix_bustypes = true) isa
          PF.ACPowerFlowData
    @test PowerFlowData(
        ACPowerFlow{NewtonRaphsonACPowerFlow}(),
        sys;
        fix_bustypes = true,
    ) isa PF.ACPowerFlowData
    @test PowerFlowData(DCPowerFlow(), sys; fix_bustypes = true) isa PF.ABAPowerFlowData
    @test PowerFlowData(PTDFDCPowerFlow(), sys; fix_bustypes = true) isa
          PF.PTDFPowerFlowData
    @test PowerFlowData(vPTDFDCPowerFlow(), sys; fix_bustypes = true) isa
          PF.vPTDFPowerFlowData
end

@testset "PowerFlowData multiperiod" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    time_steps = 24
    # TODO: "multiperiod AC still to implement"
    @test PowerFlowData(
        DCPowerFlow(),
        sys;
        time_steps = time_steps,
        fix_bustypes = true,
    ) isa PF.ABAPowerFlowData
    @test PowerFlowData(
        PTDFDCPowerFlow(),
        sys;
        time_steps = time_steps,
        fix_bustypes = true,
    ) isa
          PF.PTDFPowerFlowData
    @test PowerFlowData(
        vPTDFDCPowerFlow(),
        sys;
        time_steps = time_steps,
        fix_bustypes = true,
    ) isa
          PF.vPTDFPowerFlowData
end

@testset "System <-> PowerFlowData round trip" for ACSolver in AC_SOLVERS_TO_TEST
    # TODO currently only tested with ACPowerFlow
    # TODO test that update_system! errors if the PowerFlowData doesn't correspond to the system

    sys_original = build_system(PSISystems, "RTS_GMLC_DA_sys")
    data_original =
        PowerFlowData(ACPowerFlow{ACSolver}(), sys_original; fix_bustypes = true)

    sys_modified = deepcopy(sys_original)
    modify_rts_system!(sys_modified)
    data_modified =
        PowerFlowData(ACPowerFlow{ACSolver}(), sys_original; fix_bustypes = true)
    modify_rts_powerflow!(data_modified)

    # update_system! with unmodified PowerFlowData should result in system that yields unmodified PowerFlowData
    # (NOTE does NOT necessarily yield original system due to power redistribution)
    sys_null_updated = deepcopy(sys_original)
    PF.update_system!(sys_null_updated, data_original)
    data_null_updated =
        PowerFlowData(ACPowerFlow{ACSolver}(), sys_null_updated; fix_bustypes = true)
    @test IS.compare_values(powerflow_match_fn, data_null_updated, data_original)

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
        PowerFlowData(ACPowerFlow{ACSolver}(), sys_modified; fix_bustypes = true),
        data_modified,
    )

    # The big one: update_system! with modified PowerFlowData should result in sys_modified,
    # modulo information that is inherently lost in the PowerFlowData representation
    sys_modify_updated = deepcopy(sys_original)
    PF.update_system!(sys_modify_updated, data_modified)
    sys_mod_redist = deepcopy(sys_modified)
    PF.update_system!(
        sys_mod_redist,
        PowerFlowData(ACPowerFlow{ACSolver}(), sys_mod_redist; fix_bustypes = true),
    )
    @test IS.compare_values(powerflow_match_fn, sys_modify_updated, sys_mod_redist)
end

"""Helper function that sets availability of all sources at a given bus."""
function set_availability_at_bus(
    sys::PSY.System,
    bus::PSY.ACBus,
    availability::Bool,
)
    set_available!.(
        PSY.get_components(
            d -> !isa(d, PSY.ElectricLoad) && get_number(get_bus(d)) == get_number(bus),
            PSY.StaticInjection,
            sys,
        ),
        (availability,),
    )
end

@testset "test fix_bustypes" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14"; add_forecasts = false)
    buses = collect(PSY.get_components(PSY.ACBus, sys))
    # prevent against unknowningly having a PV bus with no available generators.
    for g in get_components(PSY.Generator, sys)
        set_available!(g, true)
    end
    first_pv = findfirst(bus -> PSY.get_bustype(bus) == PSY.ACBusTypes.PV, buses)
    pv_bus = buses[first_pv]
    set_availability_at_bus(sys, pv_bus, false)

    data = PF.PowerFlowData(PF.ACPowerFlow(), sys)
    pv_bus_row = PF.get_bus_lookup(data)[PSY.get_number(pv_bus)]
    @test PF.get_bus_type(data)[pv_bus_row, 1] == PSY.ACBusTypes.PV

    data_fixed = PF.PowerFlowData(PF.ACPowerFlow(), sys; fix_bustypes = true)
    pv_bus_row_fixed = PF.get_bus_lookup(data_fixed)[PSY.get_number(pv_bus)]
    @test PF.get_bus_type(data_fixed)[pv_bus_row_fixed, 1] == PSY.ACBusTypes.PQ

    @test PF.get_bus_type(data_fixed)[begin:end .!= pv_bus_row_fixed, 1] ==
          PF.get_bus_type(data)[begin:end .!= pv_bus_row, 1]
end
