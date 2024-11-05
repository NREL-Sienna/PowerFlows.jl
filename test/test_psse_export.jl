test_psse_export_dir = joinpath(TEST_FILES_DIR, "test_exports")
isdir(test_psse_export_dir) && rm(test_psse_export_dir; recursive = true)

function _log_assert(result, msg, comparison_name)
    result ||
        @error "Failed check: $(string(msg))$(isnothing(comparison_name) ? "" :  " ($comparison_name)")"
    return result
end
"""
If the expression is false, log an error; in any case, pass through the result of the
expression. Optionally accepts a name to include in the error log.
"""
macro log_assert(ex, comparison_name = nothing)
    return :(_log_assert($(esc(ex)), $(string(ex)), $(esc(comparison_name))))
end

"""
Compare the two dataframes by column. Specify tolerances using kwargs; tolerances default to
default_tol. If tolerance is `nothing`, skip that column. Otherwise, if the column is
floating point, compare element-wise with `isapprox(atol = tolerance)`; if not, test strict
equality element-wise. Optionally accepts a name to include in any failure logs.
"""
function compare_df_within_tolerance(
    comparison_name::String,
    df1::DataFrame,
    df2::DataFrame,
    default_tol = SYSTEM_REIMPORT_COMPARISON_TOLERANCE;
    kwargs...,
)
    result = true
    n_rows_match = (@log_assert size(df1, 1) == size(df2, 1) comparison_name)
    result &= n_rows_match
    result &= (@log_assert names(df1) == names(df2) comparison_name)
    result &= (@log_assert eltype.(eachcol(df1)) == eltype.(eachcol(df2)) comparison_name)
    n_rows_match || return result  # Can't compare the cols if number of rows doesn't match
    for (colname, my_eltype, col1, col2) in
        zip(names(df1), eltype.(eachcol(df1)), eachcol(df1), eachcol(df2))
        my_tol = (Symbol(colname) in keys(kwargs)) ? kwargs[Symbol(colname)] : default_tol
        isnothing(my_tol) && continue
        inner_result = if (my_eltype <: AbstractFloat)
            all(isapprox.(col1, col2; atol = my_tol))
        else
            all(IS.isequivalent.(col1, col2))
        end
        inner_result ||
            (@error "Mismatch on $colname$((my_eltype <: AbstractFloat) ? ", max discrepancy $(maximum(abs.(col2 - col1)))" : "") ($comparison_name)")
        result &= inner_result
    end
    return result
end

compare_df_within_tolerance(
    df1::DataFrame,
    df2::DataFrame,
    default_tol = SYSTEM_REIMPORT_COMPARISON_TOLERANCE;
    kwargs...,
) = compare_df_within_tolerance("unnamed", df1, df2, default_tol; kwargs...)

function compare_component_values(
    sys1::System,
    sys2::System;
    exclude_reactive_power = false,
)
    (PSY.get_units_base(sys1) == PSY.get_units_base(sys2)) ||
        throw(ArgumentError("Can only compare systems with the same units system"))
    # TODO rewrite to not depend on the old `_states` DataFrame-based functions
    reactive_power_tol =
        exclude_reactive_power ? nothing : SYSTEM_REIMPORT_COMPARISON_TOLERANCE
    result = true
    result &= compare_df_within_tolerance(
        "Bus_states",
        PF.Bus_states(sys1),
        PF.Bus_states(sys2);
        bus_name = nothing,
    )
    result &= compare_df_within_tolerance(
        "Line_states",
        PF.Line_states(sys1),
        PF.Line_states(sys2);
        line_name = nothing,
        active_flow = nothing,
        reactive_flow = nothing,
    )
    result &= compare_df_within_tolerance(
        "StaticLoad_states",
        PF.StaticLoad_states(sys1),
        PF.StaticLoad_states(sys2),
    )
    result &= compare_df_within_tolerance(
        "FixedAdmittance_states",
        PF.FixedAdmittance_states(sys1),
        PF.FixedAdmittance_states(sys2);
        load_name = nothing,
    )
    thermals1 = sort(PF.ThermalStandard_states(sys1))
    thermals2 = filter(
        :Generator_name => in(thermals1[!, :Generator_name]),
        sort(PF.ThermalStandard_states(sys2)),
    )
    result &= compare_df_within_tolerance(
        "ThermalStandard_states",
        thermals1,
        thermals2;
        rating = nothing,
        reactive_power = reactive_power_tol,
    )
    gens1 = sort(append!(PF.Generator_states(sys1), PF.Source_states(sys1)))
    gens2 = sort(append!(PF.Generator_states(sys2), PF.Source_states(sys2)))
    result &=
        compare_df_within_tolerance(
            "Generator_states_Source_states",
            gens1,
            gens2;
            rating = nothing,
            reactive_power = reactive_power_tol,
        )
    result &= compare_df_within_tolerance(
        "Transformer2W_states",
        PF.Transformer2W_states(sys1),
        PF.Transformer2W_states(sys2);
        Transformer_name = nothing,
        active_power_flow = nothing,
        reactive_power_flow = nothing,
    )
    result &= compare_df_within_tolerance(
        "TapTransformer_states",
        PF.TapTransformer_states(sys1),
        PF.TapTransformer_states(sys2);
        r = 1e-4, x = 0.005,
    )
    result &= compare_df_within_tolerance(
        "FixedAdmittance_states",
        PF.FixedAdmittance_states(sys1),
        PF.FixedAdmittance_states(sys2);
        load_name = nothing,
    )
    return result
end

# If we have a name like "Bus1-Bus2-OtherInfo," reverse it to "Bus2-Bus1-OtherInfo"
function reverse_composite_name(name::String)
    parts = split(name, "-")
    (length(parts) > 2) || return name
    return join([parts[2], parts[1], parts[3:end]...], "-")
end

loose_system_match_fn(a::Float64, b::Float64) =
    isapprox(a, b; atol = SYSTEM_REIMPORT_COMPARISON_TOLERANCE) || IS.isequivalent(a, b)
loose_system_match_fn(a, b) = IS.isequivalent(a, b)

function compare_systems_loosely(sys1::PSY.System, sys2::PSY.System;
    bus_name_mapping = Dict{String, String}(),
    # TODO when possible, also include: PSY.FixedAdmittance, PSY.Arc
    include_types = [
        PSY.ACBus,
        PSY.Area,
        PSY.Line,
        PSY.LoadZone,
        PSY.StandardLoad,
        PSY.Transformer2W,
        PSY.ThermalStandard,
    ],
    # TODO when possible, don't exclude so many fields
    exclude_fields = Set([
        :ext,
        :active_power_flow,
        :reactive_power_flow,
        :ramp_limits,
        :time_limits,
        :services,
        :angle_limits,
    ]),
    exclude_fields_for_type = Dict(
        PSY.ThermalStandard => Set([
            :prime_mover_type,
            :rating,
            :fuel,
            :dynamic_injector,
            :operation_cost,
        ]),
        PSY.LoadZone => Set([
            :peak_active_power,
            :peak_reactive_power,
        ]),
    ),
    ignore_name_order = true,
    ignore_extra_of_type = Union{PSY.Generator, PSY.StaticLoad},
    exclude_reactive_power = false)
    result = true
    exclude_reactive_power && push!(exclude_fields, :reactive_power)

    # Compare everything about the systems except the actual components
    result &= IS.compare_values(sys1, sys2; exclude = [:data])

    # Compare the components by concrete type
    for my_type in include_types
        !isconcretetype(my_type) &&
            throw(ArgumentError("All `include_types` must be concrete, got $my_type"))

        names1 = collect(PSY.get_name.(PSY.get_components(my_type, sys1)))
        predicted_names2 = replace.(names1, bus_name_mapping...)
        actual_names2 = collect(PSY.get_name.(PSY.get_components(my_type, sys2)))

        if ignore_name_order
            for (i, predicted) in enumerate(predicted_names2)
                if !(predicted in actual_names2) &&
                   reverse_composite_name(predicted) in actual_names2
                    @info "Reversing name $predicted"
                    predicted_names2[i] = reverse_composite_name(predicted)
                end
            end
        end

        if my_type <: ignore_extra_of_type
            if !isempty(setdiff(predicted_names2, actual_names2))
                @error "Predicting generator names that do not exist for $my_type"
                result = false
            end
            (Set(predicted_names2) != Set(actual_names2)) &&
                @warn "Predicted $my_type names are a strict subset of actual $my_type names"
        else
            if Set(predicted_names2) != Set(actual_names2)
                @error "Predicted names do not match actual names for $my_type"
                result = false
            end
        end

        my_excludes =
            union(Set(exclude_fields), get(exclude_fields_for_type, my_type, Set()))
        for (name1, name2) in zip(names1, predicted_names2)
            (name2 in actual_names2) || continue
            comp1 = PSY.get_component(my_type, sys1, name1)
            comp2 = PSY.get_component(my_type, sys2, name2)
            @assert !isnothing(comp2) comp2

            comparison = IS.compare_values(
                loose_system_match_fn,
                comp1,
                comp2;
                exclude = my_excludes,
            )
            result &= comparison
            if !comparison
                @error "Mismatched component LHS: $comp1"
                @error "Mismatched component RHS: $comp2"
            end
        end
    end
    return result
end

# We currently have two imperfect methods of comparing systems. TODO at some point combine into one good method
function compare_systems_wrapper(sys1::System, sys2::System, sys2_metadata = nothing;
    exclude_reactive_power = false)
    first_result = compare_component_values(sys1, sys2;
        exclude_reactive_power = exclude_reactive_power)
    second_result = compare_systems_loosely(sys1, sys2;
        exclude_reactive_power = exclude_reactive_power)
    return first_result && second_result
end

function test_power_flow(sys1::System, sys2::System; exclude_reactive_flow = false)
    result1 = solve_powerflow(ACPowerFlow(), sys1)
    result2 = solve_powerflow(ACPowerFlow(), sys2)
    reactive_power_tol =
        exclude_reactive_flow ? nothing : POWERFLOW_COMPARISON_TOLERANCE
    @test compare_df_within_tolerance("bus_results", result1["bus_results"],
        result2["bus_results"], POWERFLOW_COMPARISON_TOLERANCE)
    @test compare_df_within_tolerance("flow_results",
        sort(result1["flow_results"], names(result1["flow_results"])[2:end]),
        sort(result2["flow_results"], names(result2["flow_results"])[2:end]),
        POWERFLOW_COMPARISON_TOLERANCE; line_name = nothing, Q_to_from = reactive_power_tol,
        Q_from_to = reactive_power_tol, Q_losses = reactive_power_tol)
end

function read_system_and_metadata(raw_path, metadata_path)
    md = JSON3.read(metadata_path, Dict)
    sys = System(raw_path, md)
    return sys, md
end

read_system_and_metadata(export_subdir) = read_system_and_metadata(
    get_psse_export_paths(export_subdir)...)

function test_psse_round_trip(
    sys::System,
    exporter::PSSEExporter,
    scenario_name::AbstractString,
    export_location::AbstractString;
    do_power_flow_test = true,
    exclude_reactive_flow = false,
)
    raw_path, metadata_path =
        get_psse_export_paths(joinpath(export_location, scenario_name))
    @test !isfile(raw_path)
    @test !isfile(metadata_path)

    write_export(exporter, scenario_name)
    @test isfile(raw_path)
    @test isfile(metadata_path)

    sys2, sys2_metadata = read_system_and_metadata(raw_path, metadata_path)
    @test compare_systems_wrapper(sys, sys2, sys2_metadata)
    do_power_flow_test &&
        test_power_flow(sys, sys2; exclude_reactive_flow = exclude_reactive_flow)
end

"Test that the two raw files are exactly identical and the two metadata files parse to identical JSON"
function test_psse_export_strict_equality(
    raw1,
    metadata1,
    raw2,
    metadata2;
    exclude_metadata_keys = ["case_name"],
)
    open(raw1, "r") do handle1
        open(raw2, "r") do handle2
            @test countlines(handle1) == countlines(handle2)
            for (line1, line2) in zip(readlines(handle1), readlines(handle2))
                @test line1 == line2
            end
        end
    end

    parsed1 = JSON3.read(metadata1, Dict)
    parsed2 = JSON3.read(metadata2, Dict)
    for key in exclude_metadata_keys
        parsed1[key] = nothing
        parsed2[key] = nothing
    end
    @test parsed1 == parsed2
end

function load_test_system()
    # TODO commit to either providing this file or not requiring it
    data_dir = joinpath(dirname(dirname(dirname(pathof(PowerFlows)))), "pf_data")
    sys_file = joinpath(data_dir, "twofortybus", "Marenas", "system_240[32].json")
    if !isfile(sys_file)
        @warn "Skipping test with system_240[32].json, file does not exist"
        return
    end
    sys = with_logger(SimpleLogger(Error)) do
        System(sys_file)
    end
    set_units_base_system!(sys, UnitSystem.SYSTEM_BASE)
    return sys
end

# I test so much, my tests have tests
@testset "Test system comparison utilities" begin
    sys = load_test_system()
    isnothing(sys) && return

    @test compare_systems_wrapper(sys, sys)
    @test compare_systems_wrapper(sys, deepcopy(sys))
end

@testset "PSSE Exporter with system_240[32].json, v33" begin
    sys = load_test_system()
    isnothing(sys) && return

    # PSS/E version must be one of the supported ones
    @test_throws ArgumentError PSSEExporter(sys, :vNonexistent, test_psse_export_dir)

    # Reimported export should be comparable to original system
    export_location = joinpath(test_psse_export_dir, "v33", "system_240")
    exporter = PSSEExporter(sys, :v33, export_location)
    test_psse_round_trip(sys, exporter, "basic", export_location;
        exclude_reactive_flow = true)  # TODO why is reactive flow not matching?

    # Exporting the exact same thing again should result in the exact same files
    write_export(exporter, "basic2")
    test_psse_export_strict_equality(
        get_psse_export_paths(joinpath(export_location, "basic"))...,
        get_psse_export_paths(joinpath(export_location, "basic2"))...)

    # Updating with a completely different system should fail
    different_system = build_system(PSITestSystems, "c_sys5_all_components")
    @test_throws ArgumentError update_exporter!(exporter, different_system)

    # Updating with the exact same system should result in the exact same files
    update_exporter!(exporter, sys)
    write_export(exporter, "basic3")
    test_psse_export_strict_equality(
        get_psse_export_paths(joinpath(export_location, "basic"))...,
        get_psse_export_paths(joinpath(export_location, "basic3"))...)

    # Updating with changed value should result in a different reimport (System version)
    sys2 = deepcopy(sys)
    line_to_change = first(get_components(Line, sys2))
    set_rating!(line_to_change, get_rating(line_to_change) * 12345.6)
    update_exporter!(exporter, sys2)
    write_export(exporter, "basic4")
    reread_sys2, sys2_metadata =
        read_system_and_metadata(joinpath(export_location, "basic4"))
    @test compare_systems_wrapper(sys2, reread_sys2, sys2_metadata)
    @test_logs((:error, r"Mismatch on rate"), (:error, r"values do not match"),
        match_mode = :any, min_level = Logging.Error,
        compare_systems_wrapper(sys, reread_sys2, sys2_metadata))
    test_power_flow(sys2, reread_sys2; exclude_reactive_flow = true)  # TODO why is reactive flow not matching?
end

@testset "PSSE Exporter with RTS_GMLC_DA_sys, v33" begin
    sys = create_pf_friendly_rts_gmlc()
    set_units_base_system!(sys, UnitSystem.SYSTEM_BASE)

    # PSS/E version must be one of the supported ones
    @test_throws ArgumentError PSSEExporter(sys, :vNonexistent, test_psse_export_dir)

    # Reimported export should be comparable to original system
    export_location = joinpath(test_psse_export_dir, "v33", "rts_gmlc")
    exporter = PSSEExporter(sys, :v33, export_location)
    test_psse_round_trip(sys, exporter, "basic", export_location;
        exclude_reactive_flow = true)  # TODO why is reactive flow not matching?

    # Exporting the exact same thing again should result in the exact same files
    write_export(exporter, "basic2")
    test_psse_export_strict_equality(
        get_psse_export_paths(joinpath(export_location, "basic"))...,
        get_psse_export_paths(joinpath(export_location, "basic2"))...)

    # Updating with a completely different system should fail
    different_system = build_system(PSITestSystems, "c_sys5_all_components")
    @test_throws ArgumentError update_exporter!(exporter, different_system)

    # Updating with the exact same system should result in the exact same files
    update_exporter!(exporter, sys)
    write_export(exporter, "basic3")
    test_psse_export_strict_equality(
        get_psse_export_paths(joinpath(export_location, "basic"))...,
        get_psse_export_paths(joinpath(export_location, "basic3"))...)

    # Updating with changed value should result in a different reimport (System version)
    sys2 = deepcopy(sys)
    modify_rts_system!(sys2)
    update_exporter!(exporter, sys2)
    write_export(exporter, "basic4")
    reread_sys2, sys2_metadata =
        read_system_and_metadata(joinpath(export_location, "basic4"))
    @test compare_systems_wrapper(sys2, reread_sys2, sys2_metadata)
    @test_logs((:error, r"values do not match"),
        (:error, r"Mismatch on active_power"), (:error, r"Mismatch on reactive_power"),
        (:error, r"Mismatch on Vm"), (:error, r"Mismatch on θ"),
        match_mode = :any, min_level = Logging.Error,
        compare_systems_wrapper(sys, reread_sys2, sys2_metadata))
    test_power_flow(sys2, reread_sys2; exclude_reactive_flow = true)  # TODO why is reactive flow not matching?

    # Updating with changed value should result in a different reimport (PowerFlowData version)
    exporter = PSSEExporter(sys, :v33, export_location)
    pf2 = PowerFlowData(ACPowerFlow(), sys)
    # This modifies the PowerFlowData in the same way that modify_rts_system! modifies the
    # system, so the reimport should be comparable to sys2 from above
    modify_rts_powerflow!(pf2)
    update_exporter!(exporter, pf2)
    write_export(exporter, "basic5")
    reread_sys3, sys3_metadata =
        read_system_and_metadata(joinpath(export_location, "basic5"))
    @test compare_systems_wrapper(sys2, reread_sys3, sys3_metadata;
        exclude_reactive_power = true)  # TODO why is reactive power not matching?
    @test_logs((:error, r"values do not match"),
        (:error, r"Mismatch on active_power"), (:error, r"Mismatch on reactive_power"),
        (:error, r"Mismatch on Vm"), (:error, r"Mismatch on θ"),
        match_mode = :any, min_level = Logging.Error,
        compare_systems_wrapper(sys, reread_sys3, sys3_metadata))
    test_power_flow(sys2, reread_sys3; exclude_reactive_flow = true)  # TODO why is reactive flow not matching?

    # Exporting with write_comments should be comparable to original system
    exporter = PSSEExporter(sys, :v33, export_location; write_comments = true)
    test_psse_round_trip(sys, exporter, "basic6", export_location;
        exclude_reactive_flow = true)  # TODO why is reactive flow not matching?
end

@testset "Test exporter helper functions" begin
    @test PF._psse_bus_numbers([2, 3, 999_997, 999_998, 1_000_001, 1]) ==
          Dict(
        2 => 2,
        3 => 3,
        999_997 => 999_997,
        999_998 => 899_998,
        1_000_001 => 4,
        1 => 1,
    )
    @test !PF._is_valid_psse_name("a pretty long name")
    @test !PF._is_valid_psse_name("-bad")
    @test PF._is_valid_psse_name(raw"¯\_(ツ)_/¯")
    @test PF._psse_bus_names(["-bad1", "compliant", "BUS_100", "-bad2", "ok just too long"],
        [10, 2, 3, 4, 5], Dict(10 => 100, 2 => 20, 3 => 30, 4 => 40, 5 => 50)) ==
          Dict("-bad1" => "BUS_100-", "compliant" => "compliant", "BUS_100" => "BUS_100",
        "-bad2" => "BUS_40", "ok just too long" => "ok just too ")
    @test PF.create_component_ids(
        ["generator-1234-AB", "123_CT_7", "load1234", "load1334"], [1, 1, 2, 2]) ==
          Dict((1, "generator-1234-AB") => "AB", (1, "123_CT_7") => "7",
        (2, "load1234") => "34", (2, "load1334") => "35")
end

# TODO add tests for unit system agnosticism
