test_psse_export_dir = joinpath(TEST_FILES_DIR, "test_psse_exports")  # at some point could move this to temp files
isdir(test_psse_export_dir) && rm(test_psse_export_dir; recursive = true)

"""
Compare the two dataframes by column. Specify tolerances using kwargs; tolerances default to
default_tol. If tolerance is `nothing`, skip that column. Otherwise, if the column is
floating point, compare element-wise with `isapprox(atol = tolerance)`; if not, test strict
equality element-wise.
"""
function test_diff_within_tolerance(
    df1::DataFrame,
    df2::DataFrame,
    default_tol = PF.SYSTEM_EXPORT_TOLERANCE;
    kwargs...,
)
    @test names(df1) == names(df2)
    @test eltype.(eachcol(df1)) == eltype.(eachcol(df2))
    for (colname, my_eltype, col1, col2) in
        zip(names(df1), eltype.(eachcol(df1)), eachcol(df1), eachcol(df2))
        my_tol = (Symbol(colname) in keys(kwargs)) ? kwargs[Symbol(colname)] : default_tol
        isnothing(my_tol) && continue
        success = if my_eltype <: AbstractFloat
            @test all(isapprox.(col1, col2; atol = my_tol))
        else
            @test all(isequal.(col1, col2))
        end
        (success isa Test.Pass) || @error "Mismatch on $colname"
    end
end

function compare_component_values(sys1::System, sys2::System)
    # TODO rewrite to not depend on the old `_states` DataFrame-based functions
    test_diff_within_tolerance(PF.Bus_states(sys1), PF.Bus_states(sys2); bus_name = nothing)
    test_diff_within_tolerance(
        PF.Line_states(sys1),
        PF.Line_states(sys2);
        line_name = nothing,
        active_flow = nothing,
        reactive_flow = nothing,
    )
    test_diff_within_tolerance(PF.StandardLoad_states(sys1), PF.StandardLoad_states(sys2))
    test_diff_within_tolerance(
        PF.FixedAdmittance_states(sys1),
        PF.FixedAdmittance_states(sys2);
        load_name = nothing,
    )
    thermals1 = sort(PF.ThermalStandard_states(sys1))
    thermals2 = filter(
        :Generator_name => in(thermals1[!, :Generator_name]),
        sort(PF.ThermalStandard_states(sys2)),
    )
    test_diff_within_tolerance(thermals1, thermals2; rating = nothing)
    gens1 = sort(append!(PF.Generator_states(sys1), PF.Source_states(sys1)))
    gens2 = sort(PF.Generator_states(sys2))
    test_diff_within_tolerance(gens1, gens2; rating = nothing, Generator_name = nothing)
    test_diff_within_tolerance(
        PF.Transformer2W_states(sys1),
        PF.Transformer2W_states(sys2);
        Transformer_name = nothing,
        active_power_flow = nothing,
        reactive_power_flow = nothing,
    )
    test_diff_within_tolerance(
        PF.TapTransformer_states(sys1),
        PF.TapTransformer_states(sys2),
    )
    test_diff_within_tolerance(
        PF.FixedAdmittance_states(sys1),
        PF.FixedAdmittance_states(sys2);
        load_name = nothing,
    )
end

# If we have a name like "Bus1-Bus2-OtherInfo," reverse it to "Bus2-Bus1-OtherInfo"
function reverse_composite_name(name::String)
    parts = split(name, "-")
    (length(parts) > 2) || return name
    return join([parts[2], parts[1], parts[3:end]...], "-")
end

loose_system_match_fn(a::Float64, b::Float64) =
    isapprox(a, b; atol = PF.SYSTEM_EXPORT_TOLERANCE) || IS.isequivalent(a, b)
loose_system_match_fn(a, b) = IS.isequivalent(a, b)

function compare_systems_loosely(sys1::PSY.System, sys2::PSY.System;
    bus_name_mapping = Dict{String, String}(),
    # TODO when possible, also include: PSY.FixedAdmittance
    include_types = [
        PSY.ACBus,
        PSY.Arc,
        PSY.Area,
        PSY.Line,
        PSY.LoadZone,
        PSY.StandardLoad,
        PSY.Transformer2W,
        PSY.ThermalStandard,
    ],
    # TODO when possible, don't exclude: :bustype, :number, :angle, :magnitude, :variable, maybe more
    exclude_fields = Set([
        :name,
        :ext,
        :bustype,
        :angle,
        :magnitude,
        :active_power_flow,
        :reactive_power_flow,
        :internal,
    ]),
    # TODO when possible, don't exclude these things
    exclude_fields_for_type = Dict(
        PSY.ThermalStandard => Set([
            :prime_mover_type,
            :rating,
            :fuel,
            :active_power_limits,
            :reactive_power_limits,
            :dynamic_injector,
            :operation_cost,
        ]),
    ),
    ignore_name_order = true,
    ignore_extra_gens = true)
    result = true

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

        if ignore_extra_gens && my_type <: PSY.Generator
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
                comp1,
                comp2;
                exclude = my_excludes,
                match_fn = loose_system_match_fn,
            )
            result &= comparison
            if !comparison
                @show comp1
                @show comp2
            end
        end
    end
    return result
end

# We currently have two imperfect methods of comparing systems. TODO at some point combine into one good method
function compare_systems_wrapper(sys1::System, sys2::System, sys2_metadata)
    compare_component_values(sys1, sys2)
    compare_systems_loosely(
        sys1,
        sys2;
        bus_name_mapping = sys2_metadata["Bus_Name_Mapping"],
    )
end

function test_psse_round_trip(
    sys::System,
    exporter::PSSEExporter,
    scenario_name::AbstractString,
    year::Int,
    export_location::AbstractString,
)
    raw_path, metadata_path = get_paths(scenario_name, year, export_location)
    @test !isfile(raw_path)
    @test !isfile(metadata_path)

    write_export(exporter, scenario_name, year, export_location)
    @test isfile(raw_path)
    @test isfile(metadata_path)

    sys2 = System(raw_path)
    sys2_metadata = PF.JSON.parsefile(metadata_path)

    set_units_base_system!(sys, UnitSystem.SYSTEM_BASE)
    set_units_base_system!(sys2, UnitSystem.SYSTEM_BASE)

    compare_systems_wrapper(sys, sys2, sys2_metadata)
end

@testset "PSSE Exporter with system_240[32].json, v33" begin
    # TODO commit to either providing this file or not requiring it
    sys_file = joinpath(PF.DATA_DIR, "twofortybus", "Marenas", "system_240[32].json")
    if !isfile(sys_file)
        @warn "Skipping test with system_240[32].json, file does not exist"
        return
    end
    sys = with_logger(SimpleLogger(Error)) do
        System(sys_file)
    end
    @test_throws ArgumentError PSSEExporter(sys, :vNonexistent)

    exporter_33 = PSSEExporter(sys, :v33)
    export_location = joinpath(test_psse_export_dir, "v33", "system_240")
    test_psse_round_trip(sys, exporter_33, "basic", 2024, export_location)
end

# TODO make this pass
# @testset "PSSE Exporter with RTS_GMLC_DA_sys, v33" begin
#     sys = build_system(PSISystems, "RTS_GMLC_DA_sys")
#     @test_throws ArgumentError PSSEExporter(sys, :vNonexistent)

#     exporter_33 = PSSEExporter(sys, :v33)
#     export_location = joinpath(test_psse_export_dir, "v33", "RTS_GMLC_DA_sys")
#     test_psse_round_trip(sys, exporter_33, "basic", 2024, export_location)
# end

# TODO v44
