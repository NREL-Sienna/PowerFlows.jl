const PSSE_EXPORT_SUPPORTED_VERSIONS = [:v33]  # TODO add :v34
const PSSE_DEFAULT = ""  # Used below in cases where we want to insert an empty field to signify the PSSE default
const PSSE_BUS_TYPE_MAP = Dict(
    PSY.ACBusTypes.PQ => 1,
    PSY.ACBusTypes.PV => 2,
    PSY.ACBusTypes.REF => 3,
    PSY.ACBusTypes.SLACK => 3,
    PSY.ACBusTypes.ISOLATED => 4,
)
# Splat this out where the eight ownership fields are necessary
const PSSE_DEFAULT_OWNERSHIP = let
    O1, O2, O3, O4 = PSSE_DEFAULT, PSSE_DEFAULT, PSSE_DEFAULT, PSSE_DEFAULT
    F1, F2, F3, F4 = PSSE_DEFAULT, PSSE_DEFAULT, PSSE_DEFAULT, PSSE_DEFAULT
    O1, F1, O2, F2, O3, F3, O4, F4
end

# TODO consider adding this to IS
"""
A context manager similar to `Logging.with_logger` that sets the system's units to the given
value, executes the function, then sets them back. Suppresses logging below `Warn` from
internal calls to `set_units_base_system!`.
"""
function with_units(f::Function, sys::System, units::Union{PSY.UnitSystem, String})
    old_units = PSY.get_units_base(sys)
    Logging.with_logger(Logging.SimpleLogger(Logging.Warn)) do
        PSY.set_units_base_system!(sys, units)
    end
    try
        f()
    finally
        Logging.with_logger(Logging.SimpleLogger(Logging.Warn)) do
            PSY.set_units_base_system!(sys, old_units)
        end
    end
end

"""
Structure to perform an export from a Sienna System, plus optional updates from
`PowerFlowData`, to the PSS/E format. Construct from a `System` and a PSS/E version, update
using `update_exporter` with any new data as relevant, and perform the export with
`write_export`.

# Arguments:
  - `base_system::PSY.System`: the system to be exported. Later updates may change power
    flow-related values but may not fundamentally alter the system
  - `psse_version::Symbol`: the version of PSS/E to target, must be one of
    `PSSE_EXPORT_SUPPORTED_VERSIONS`
"""
mutable struct PSSEExporter
    # Internal fields are very much subject to change as I iterate on the best way to do
    # this! For instance, the final version will almost certainly not store an entire System
    system::PSY.System
    psse_version::Symbol

    function PSSEExporter(base_system::PSY.System, psse_version::Symbol)
        (psse_version in PSSE_EXPORT_SUPPORTED_VERSIONS) ||
            throw(
                ArgumentError(
                    "PSS/E version $psse_version is not supported, must be one of $PSSE_EXPORT_SUPPORTED_VERSIONS",
                ),
            )
        system = deepcopy(base_system)
        new(system, psse_version)
    end
end

function _validate_same_system(sys1::PSY.System, sys2::PSY.System)
    return IS.get_uuid(PSY.get_internal(sys1)) == IS.get_uuid(PSY.get_internal(sys2))
end

"""
Update the `PSSEExporter` with new `data`.

# Arguments:
  - `exporter::PSSEExporter`: the exporter to update
  - `data::PSY.PowerFlowData`: the new data. Must correspond to the `System` with which the
    exporter was constructor
"""
function update_exporter!(exporter::PSSEExporter, data::PowerFlowData)
    # NOTE this relies on exporter.system being a deepcopy of the original system so we're not changing that one here
    update_system!(exporter.system, data)
end

# TODO solidify the notion of sameness we care about here
"""
Update the `PSSEExporter` with new `data`.

# Arguments:
  - `exporter::PSSEExporter`: the exporter to update
  - `data::PSY.System`: system containing the new data. Must be fundamentally the same
  `System` as the one with which the exporter was constructed, just with different values
"""
function update_exporter!(exporter::PSSEExporter, data::PSY.System)
    _validate_same_system(exporter.system, data) || throw(
        ArgumentError(
            "System passed to update_exporter must be the same system as the one with which the exporter was constructed, just with different values",
        ),
    )
    exporter.system = deepcopy(data)
end

"""
`join` with a newline at the end, delimeter defaults to \", \". If `strip_trailing_empties`,
trailing entries of `iterator` equal to the empty string are ignored.
"""
function joinln(io::IO, iterator, delim = ", "; strip_trailing_empties = true)  # TODO maybe remove the space from the delim?
    (strip_trailing_empties && ("" in iterator)) &&
        (iterator = iterator[1:findlast(!=(""), iterator)])
    join(io, iterator, delim)
    println(io)
end

_permissive_parse_int(x) = Int64(parse(Float64, x))  # Parses "1.0" as 1, errors on "1.5"

_psse_quote_string(s::String) = "'$s'"

branch_to_bus_numbers(branch) =
    (PSY.get_number.((PSY.get_from_bus(branch), PSY.get_to_bus(branch))))

"Throw a `NotImplementedError` if the `psse_version` is not `:v33`"
check_33(psse_version::Symbol) =
    (psse_version == :v33) ||
    throw(IS.NotImplementedError("Only implemented for psse_version $(:v33)"))

function _validate_container_number(unparsed::String)
    parsed = try
        _permissive_parse_int(unparsed)
    catch e
        if e isa Union{InexactError, ArgumentError}
            throw(ArgumentError("container name $unparsed could not be parsed as an integer"))
        else
            rethrow(e)
        end
    end
    (parsed in 1:9999) && (return parsed)
    throw(ArgumentError("container number $parsed is out of range"))
end

"Validate that the Sienna area/zone names parse as PSS/E-compatible area/zone numbers, output a mapping"
_psse_container_numbers(container_names::Vector{String}) =
    DataStructures.OrderedDict(
        name => _validate_container_number(name) for name in container_names
    )

"WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Case Identification Data"
function _write_case_identification_data(
    io::IO,
    md::AbstractDict,
    sys::System,
    psse_version::Symbol,
    case_name::String,
)
    check_33(psse_version)

    # Record 1
    IC = 0
    SBASE = PSY.get_base_power(sys)
    REV = 33
    XFRRAT = 0
    NXFRAT = 1  # TODO why?
    BASFRQ = PSY.get_frequency(sys)
    joinln(io, [IC, SBASE, REV, XFRRAT, NXFRAT, BASFRQ])

    # Record 2
    (length(case_name) <= 60) ||
        throw(ArgumentError("case_name may be up to 60 characters"))
    println(io, case_name)

    # Record 3
    line3 = "Created by PowerFlows.jl on $(Dates.now())"
    @assert length(line3) <= 60
    println(io, line3)
    md["record_groups"]["Case Identification Data"] = true
end

"""
Given a vector of Sienna bus numbers, create a dictionary from Sienna bus number to
PSS/E-compatible bus number. Assumes that the Sienna bus numbers are positive and unique.
Guarantees determinism: if the input contains the same bus numbers in the same order, the
output will. Guarantees minimal changes: that if an existing bus number is compliant, it
will not be changed.

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Bus Data
"""
function _psse_bus_numbers(buses::Vector{Int64})
    used_numbers = Set{Int64}()
    sizehint!(used_numbers, length(buses))
    mapping = OrderedDict{Int64, Int64}()
    sizehint!(mapping, length(buses))

    for original_number in buses
        if original_number in 1:999_997
            mapping[original_number] = original_number
            push!(used_numbers, original_number)
        end
    end
    for original_number in buses
        haskey(mapping, original_number) && continue
        new_number = original_number
        new_number %= 1_000_000
        (new_number in 999_997:1_000_000) && (new_number -= 100_000)
        while new_number in used_numbers
            new_number += 1
        end
        mapping[original_number] = new_number
        push!(used_numbers, new_number)
    end
    return mapping
end

function _is_valid_psse_name(name::String)
    (length(name) <= 12) || (return false)
    (length(name) >= 1) && (first(name) == '-') && (return false)
    return true  # Does the allowance for special characters cover *any* special characters?
end

"""
Given a vector of Sienna bus names, create a dictionary from Sienna bus name to
PSS/E-compatible bus name. Guarantees determinism and minimal changes.

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Bus Data
"""
function _psse_bus_names(
    buses::Vector{String},
    bus_numbers::Vector{Int64},
    bus_number_mapping::AbstractDict{Int64, Int64},
)
    used_names = Set{String}()
    sizehint!(used_names, length(buses))
    mapping = OrderedDict{String, String}()
    sizehint!(mapping, length(buses))

    for original_name in buses
        if _is_valid_psse_name(original_name)
            mapping[original_name] = original_name
            push!(used_names, original_name)
        end
    end
    for (original_name, original_number) in zip(buses, bus_numbers)
        haskey(mapping, original_name) && continue
        new_name = first(original_name, 12)
        if !_is_valid_psse_name(new_name) || new_name in used_names
            new_name = "BUS_$(bus_number_mapping[original_number])"
            while new_name in used_names
                new_name *= "-"
            end
        end
        @assert _is_valid_psse_name(new_name) new_name
        mapping[original_name] = new_name
        push!(used_names, new_name)
    end
    return mapping
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Bus Data. Sienna voltage limits treated as PSS/E
normal voltage limits; PSSE emergency voltage limits left as default.
"""
function _write_bus_data(io::IO, md::AbstractDict, sys::System, psse_version::Symbol)
    check_33(psse_version)

    buses = sort!(collect(PSY.get_components(PSY.Bus, sys)); by = PSY.get_number)
    old_bus_numbers = PSY.get_number.(buses)
    bus_number_mapping = _psse_bus_numbers(old_bus_numbers)
    bus_name_mapping =
        _psse_bus_names(PSY.get_name.(buses), old_bus_numbers, bus_number_mapping)
    for bus in buses
        I = bus_number_mapping[PSY.get_number(bus)]
        NAME = _psse_quote_string(bus_name_mapping[PSY.get_name(bus)])
        BASKV = PSY.get_base_voltage(bus)
        IDE = PSSE_BUS_TYPE_MAP[PSY.get_bustype(bus)]
        AREA = md["area_mapping"][PSY.get_name(PSY.get_area(bus))]
        ZONE = md["zone_number_mapping"][PSY.get_name(PSY.get_load_zone(bus))]
        OWNER = PSSE_DEFAULT
        VM = PSY.get_magnitude(bus)
        VA = rad2deg(PSY.get_angle(bus))
        NVHI = PSY.get_voltage_limits(bus).max
        NVLO = PSY.get_voltage_limits(bus).min
        EVHI = PSSE_DEFAULT
        EVLO = PSSE_DEFAULT
        joinln(io, [I, NAME, BASKV, IDE, AREA, ZONE, OWNER, VM, VA, NVHI, NVLO, EVHI, EVLO])
    end
    println(io, "0")  # End of section

    md["bus_number_mapping"] = bus_number_mapping
    md["bus_name_mapping"] = bus_name_mapping
    md["record_groups"]["Bus Data"] = true
end

function _increment_component_char(component_char::Char)
    (component_char == '9') && return 'A'
    (component_char == 'Z') && return '0'
    return component_char + 1
end

function _increment_component_id(component_id::String)
    carry = (last(component_id) == 'Z')
    if length(component_id) == 1
        carry && return '0' * _increment_component_char(first(component_id))
        return string(_increment_component_char(first(component_id)))
    end
    return (carry ? _increment_component_char(first(component_id)) : first(component_id)) *
           _increment_component_char(last(component_id))
end

"""
Try to make an informative one or two character name for the load/generator/etc.

  - "generator-1234-AB" -> "AB"
  - "123_CT_7" -> "7"
  - "load1234" -> "34"
"""
function _first_choice_gen_id(name::String)
    my_split = argmax(length, [split(name, "-"), split(name, "_")])
    return uppercase(last(last(my_split), 2))
end

"""
Given a vector of component names and a corresponding vector of container IDs (e.g., bus
numbers), create unique-per-container PSS/E-compatible IDs, output a dictionary from
(container ID, component name) to PSS/E-compatible component ID. The "singles_to_1" flag
detects components that are the only one on their bus and gives them the name "1".
"""
function create_component_ids(
    component_names::Vector{<:String},
    container_ids::Vector{T};
    singles_to_1 = false,
) where {T}
    id_mapping = Dict{Tuple{T, String}, String}()
    sizehint!(id_mapping, length(component_names))
    ids_by_container = Dict{T, Vector{String}}()

    for (name, container_id) in zip(component_names, container_ids)
        haskey(ids_by_container, container_id) ||
            (ids_by_container[container_id] = Vector{String}())
        my_blocked = ids_by_container[container_id]
        id = _first_choice_gen_id(name)
        while id in my_blocked
            id = _increment_component_id(id)
        end
        id_mapping[(container_id, name)] = id
        push!(my_blocked, id)
    end
    if singles_to_1
        for (name, container_id) in zip(component_names, container_ids)
            (length(ids_by_container[container_id]) == 1) &&
                (id_mapping[(container_id, name)] = "1")
        end
    end
    return id_mapping
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Load Data
"""
function _write_load_data(io::IO, md::AbstractDict, sys::System, psse_version::Symbol)
    check_33(psse_version)

    loads = sort!(collect(PSY.get_components(PSY.StaticLoad, sys)); by = PSY.get_name)
    load_name_mapping =
        create_component_ids(
            PSY.get_name.(loads),
            PSY.get_number.(PSY.get_bus.(loads));
            singles_to_1 = true,
        )
    for load in loads
        sienna_bus_number = PSY.get_number(PSY.get_bus(load))
        I = md["bus_number_mapping"][sienna_bus_number]
        ID = _psse_quote_string(load_name_mapping[(sienna_bus_number, PSY.get_name(load))])  # TODO should this be quoted?
        STATUS = PSY.get_available(load) ? 1 : 0
        AREA = PSSE_DEFAULT  # defaults to bus's area
        ZONE = PSSE_DEFAULT  # defaults to zone's area
        PL = PSY.get_constant_active_power(load)
        QL = PSY.get_constant_reactive_power(load)
        IP, IQ, YP, YQ = with_units(sys, PSY.UnitSystem.DEVICE_BASE) do
            PSY.get_current_active_power(load),
            PSY.get_current_reactive_power(load),
            PSY.get_impedance_active_power(load),
            PSY.get_impedance_reactive_power(load)
        end
        OWNER = PSSE_DEFAULT  # defaults to bus's owner
        SCALE = PSSE_DEFAULT  # TODO reconsider
        INTRPT = PSSE_DEFAULT  # TODO reconsider
        joinln(
            io,
            [I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER, SCALE, INTRPT],
        )
    end
    println(io, "0")  # End of section
    md["load_name_mapping"] = load_name_mapping  # TODO reshape to be better for import
    md["record_groups"]["Load Data"] = true
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Fixed Bus Shunt Data
"""
function _write_fixed_bus_shunt_data(
    io::IO,
    md::AbstractDict,
    sys::System,
    psse_version::Symbol,
)
    check_33(psse_version)

    shunts = sort!(collect(PSY.get_components(PSY.FixedAdmittance, sys)); by = PSY.get_name)
    shunt_name_mapping =
        create_component_ids(
            PSY.get_name.(shunts),
            PSY.get_number.(PSY.get_bus.(shunts));
            singles_to_1 = true,
        )
    for shunt in shunts
        sienna_bus_number = PSY.get_number(PSY.get_bus(shunt))
        I = md["bus_number_mapping"][sienna_bus_number]
        ID =
            _psse_quote_string(shunt_name_mapping[(sienna_bus_number, PSY.get_name(shunt))])  # TODO should this be quoted?
        STATUS = PSY.get_available(shunt) ? 1 : 0
        GL = real(PSY.get_Y(shunt)) * PSY.get_base_power(sys)
        BL = imag(PSY.get_Y(shunt)) * PSY.get_base_power(sys)
        joinln(
            io,
            [I, ID, STATUS, GL, BL],
        )
    end
    println(io, "0")  # End of section
    md["shunt_name_mapping"] = shunt_name_mapping  # TODO reshape to be better for import
    md["record_groups"]["Fixed Bus Shunt Data"] = true
end

"""
If the flag `sources_as_generators` is set, export `PSY.Source` instances as PSS/E
generators in addition to `PSY.Generator`s

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Fixed Bus Shunt Data
"""
function _write_generator_data(
    io::IO,
    md::AbstractDict,
    sys::System,
    psse_version::Symbol;
    sources_as_generators = false,
)
    check_33(psse_version)

    generators::Vector{PSY.StaticInjection} =
        sort!(collect(PSY.get_components(PSY.Generator, sys)); by = PSY.get_name)
    sources_as_generators && append!(generators,
        sort!(collect(PSY.get_components(PSY.Source, sys)); by = PSY.get_name))
    generator_name_mapping =
        create_component_ids(
            PSY.get_name.(generators),
            PSY.get_number.(PSY.get_bus.(generators));
            singles_to_1 = false,
        )
    for generator in generators
        sienna_bus_number = PSY.get_number(PSY.get_bus(generator))
        I = md["bus_number_mapping"][sienna_bus_number]
        ID =
            _psse_quote_string(
                generator_name_mapping[(sienna_bus_number, PSY.get_name(generator))],
            )  # TODO should this be quoted?
        PG, QG = with_units(sys, PSY.UnitSystem.SYSTEM_BASE) do
            # Doing the conversion myself due to https://github.com/NREL-Sienna/PowerSystems.jl/issues/1164
            PSY.get_active_power(generator) * PSY.get_base_power(sys),
            PSY.get_reactive_power(generator) * PSY.get_base_power(sys)
        end  # TODO fix units
        QT = PSY.get_reactive_power_limits(generator).max
        isfinite(QT) || (QT = PSSE_DEFAULT)  # Catch Inf, etc.
        QB = PSY.get_reactive_power_limits(generator).min
        isfinite(QB) || (QB = PSSE_DEFAULT)
        VS = PSY.get_magnitude(PSY.get_bus(generator))  # TODO is this correct? Should this be `get_internal_voltage` for `PSY.Source`?
        IREG = get(PSY.get_ext(generator), "IREG", PSSE_DEFAULT)
        MBASE = PSY.get_base_power(generator)
        ZR, ZX = PSSE_DEFAULT, PSSE_DEFAULT
        RT, XT = PSSE_DEFAULT, PSSE_DEFAULT  # TODO?
        GTAP = PSSE_DEFAULT
        STAT = PSY.get_available(generator) ? 1 : 0
        RMPCT = PSSE_DEFAULT
        # TODO maybe have a better default here
        PT = try
            PSY.get_active_power_limits(generator).max
        catch
            PSSE_DEFAULT
        end
        PB = try
            PSY.get_active_power_limits(generator).min
        catch
            PSSE_DEFAULT
        end
        WMOD = get(PSY.get_ext(generator), "WMOD", PSSE_DEFAULT)
        WPF = get(PSY.get_ext(generator), "WPF", PSSE_DEFAULT)
        joinln(
            io,
            [I, ID, PG, QG, QT, QB, VS, IREG, MBASE, ZR, ZX, RT, XT, GTAP, STAT,
                RMPCT, PT, PB, PSSE_DEFAULT_OWNERSHIP..., WMOD, WPF],
        )
    end
    println(io, "0")  # End of section
    md["generator_name_mapping"] = generator_name_mapping  # TODO reshape to be better for import
    md["record_groups"]["Generator Data"] = true
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Fixed Bus Shunt Data
"""
function _write_non_transformer_branch_data(
    io::IO,
    md::AbstractDict,
    sys::System,
    psse_version::Symbol,
)
    check_33(psse_version)

    # TODO can/should we be more general than `Line`?
    branches = sort!(collect(PSY.get_components(PSY.Line, sys)); by = branch_to_bus_numbers)
    branch_name_mapping =
        create_component_ids(
            PSY.get_name.(branches),
            branch_to_bus_numbers.(branches);
            singles_to_1 = false,
        )

    for branch in branches
        from_n, to_n = branch_to_bus_numbers(branch)
        I = md["bus_number_mapping"][from_n]
        J = md["bus_number_mapping"][to_n]
        CKT = branch_name_mapping[((from_n, to_n), PSY.get_name(branch))]
        @assert !(first(CKT) in ['&', '@', '*'])  # Characters with a special meaning in this context
        CKT = _psse_quote_string(CKT)
        R = PSY.get_r(branch)
        X = PSY.get_x(branch)
        B = 0.0  # TODO iron out the details of B vs. BI, BJ
        RATEA =
            RATEB =
                RATEC =
                    with_units(
                        () -> PSY.get_rating(branch),
                        sys,
                        PSY.UnitSystem.NATURAL_UNITS,
                    )
        GI, BI = 0.0, PSY.get_b(branch).from
        GJ, BJ = 0.0, PSY.get_b(branch).to
        ST = PSY.get_available(branch) ? 1 : 0
        MET = PSSE_DEFAULT
        LEN = PSSE_DEFAULT

        joinln(
            io,
            [I, J, CKT, R, X, B, RATEA, RATEB, RATEC, GI, BI, GJ, BJ, ST, MET, LEN,
                PSSE_DEFAULT_OWNERSHIP...],
        )
    end
    println(io, "0")  # End of section
    md["branch_name_mapping"] = branch_name_mapping
    md["record_groups"]["Non-Transformer Branch Data"] = true
end

"""
Given a vector of Sienna transformer names, create a dictionary from Sienna transformer name
to PSS/E-compatible transformer name. Guarantees determinism and minimal changes.

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Transformer Data
"""
function _psse_transformer_names(
    transformers::Vector{String},
    bus_numbers::Vector{Tuple{Int64, Int64}},
    bus_number_mapping::AbstractDict{Int64, Int64},
    transformer_ckt_mapping,
)
    used_names = Set{String}()
    sizehint!(used_names, length(transformers))
    mapping = OrderedDict{String, String}()
    sizehint!(mapping, length(transformers))

    for original_name in transformers
        if _is_valid_psse_name(original_name)
            mapping[original_name] = original_name
            push!(used_names, original_name)
        end
    end
    for (original_name, (orig_from, orig_to)) in zip(transformers, bus_numbers)
        haskey(mapping, original_name) && continue
        ckt = transformer_ckt_mapping[((orig_from, orig_to), original_name)]
        new_name = "B$(bus_number_mapping[orig_from])-$(bus_number_mapping[orig_to])_$ckt"
        while new_name in used_names
            new_name *= "-"
        end
        @assert _is_valid_psse_name(new_name) new_name
        mapping[original_name] = new_name
        push!(used_names, new_name)
    end
    return mapping
end

# TODO support three-winding transformers
"""
Currently only supports two-winding transformers

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Fixed Bus Shunt Data
"""
function _write_transformer_data(
    io::IO,
    md::AbstractDict,
    sys::System,
    psse_version::Symbol,
)
    check_33(psse_version)
    transformer_types =
        Union{PSY.Transformer2W, PSY.TapTransformer, PSY.PhaseShiftingTransformer}
    transformers = sort!(
        collect(PSY.get_components(transformer_types, sys));
        by = branch_to_bus_numbers,
    )
    transformer_ckt_mapping =
        create_component_ids(
            PSY.get_name.(transformers),
            branch_to_bus_numbers.(transformers);
            singles_to_1 = false,
        )
    transformer_name_mapping = _psse_transformer_names(
        PSY.get_name.(transformers),
        branch_to_bus_numbers.(transformers),
        md["bus_number_mapping"],
        transformer_ckt_mapping,
    )
    for transformer in transformers
        from_n, to_n = branch_to_bus_numbers(transformer)
        I = md["bus_number_mapping"][from_n]
        J = md["bus_number_mapping"][to_n]
        K = 0  # no third winding
        CKT = transformer_ckt_mapping[((from_n, to_n), PSY.get_name(transformer))]
        @assert !(first(CKT) in ['&', '@', '*'])  # Characters with a special meaning in this context
        CKT = _psse_quote_string(CKT)
        CW = 1  # TODO
        CZ = 1  # TODO
        CM = 1  # TODO
        MAG1 = 0.0  # TODO
        MAG2 = 0.0  # TODO
        NMETR = PSSE_DEFAULT
        NAME = _psse_quote_string(transformer_name_mapping[PSY.get_name(transformer)])
        STAT = PSY.get_available(transformer) ? 1 : 0
        VECGRP = PSSE_DEFAULT

        R1_2 = PSY.get_r(transformer)
        X1_2 = PSY.get_x(transformer)
        SBASE1_2 = PSY.get_base_power(transformer)

        # TODO verify correctness
        WINDV1 = (transformer isa PSY.TapTransformer) ? PSY.get_tap(transformer) : 1.0
        NOMV1 = 0.0  # special case: identical to bus voltage
        ANG1 = if (transformer isa PSY.PhaseShiftingTransformer)
            rad2deg(PSY.get_α(transformer))
        else
            0.0
        end
        RATA1 =
            RATB1 =
                RATC1 =
                    with_units(
                        () -> PSY.get_rating(transformer),
                        sys,
                        PSY.UnitSystem.NATURAL_UNITS,
                    )
        COD1 = PSSE_DEFAULT
        CONT1 = PSSE_DEFAULT
        RMA1 = RMI1 = VMA1 = VMI1 = PSSE_DEFAULT
        NTP1 = PSSE_DEFAULT
        TAB1 = PSSE_DEFAULT
        CR1 = CX1 = PSSE_DEFAULT
        CNXA1 = PSSE_DEFAULT

        WINDV2 = 1.0
        NOMV2 = 0.0  # special case: identical to bus voltage

        joinln(
            io,
            [I, J, K, CKT, CW, CZ, CM, MAG1, MAG2, NMETR, NAME, STAT,
                PSSE_DEFAULT_OWNERSHIP..., VECGRP],
        )
        joinln(io, [R1_2, X1_2, SBASE1_2])
        joinln(
            io,
            [WINDV1, NOMV1, ANG1, RATA1, RATB1, RATC1, COD1, CONT1, RMA1, RMI1,
                VMA1, VMI1, NTP1, TAB1, CR1, CX1, CNXA1],
        )
        joinln(io, [WINDV2, NOMV2])
    end
    println(io, "0")  # End of section
    md["transformer_ckt_mapping"] = transformer_ckt_mapping
    md["transformer_name_mapping"] = transformer_name_mapping
    md["record_groups"]["Transformer Data"] = true
end

# TODO this assumption might not be valid
"""
Assumes that the Sienna zone names are already PSS/E compatible

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Zone Data
"""
function _write_zone_data(io::IO, md::AbstractDict, sys::System, psse_version::Symbol)
    check_33(psse_version)
    zone_number_mapping = md["zone_number_mapping"]
    zones = sort!(
        collect(PSY.get_components(PSY.LoadZone, sys));
        by = x -> zone_number_mapping[PSY.get_name(x)],
    )
    for zone in zones
        name = PSY.get_name(zone)
        I = zone_number_mapping[name]
        @assert _is_valid_psse_name(name) name
        ZONAME = _psse_quote_string(name)

        joinln(io, [I, ZONAME])
    end
    println(io, "0")  # End of section
    md["record_groups"]["Zone Data"] = true
end

function _write_q_record(io::IO, md::AbstractDict, sys::System, psse_version::Symbol)
    check_33(psse_version)
    println(io, "Q")  # End of file
    md["record_groups"]["Q Record"] = true
end

function _write_skip_group(
    io::IO,
    md::AbstractDict,
    sys::System,
    psse_version::Symbol,
    section_name::String,
)
    check_33(psse_version)
    println(io, "0")
    haskey(md, "skipped_groups") || (md["skipped_groups"] = Vector{String}())
    md["record_groups"][section_name] = false
end

"Peform an export from the data contained in a `PSSEExporter` to the PSS/E file format."
function write_export(
    exporter::PSSEExporter,
    scenario_name::AbstractString,
    year::Int,
    export_location::AbstractString,
)
    # Construct paths
    export_dir = joinpath(export_location, "Raw_Export", scenario_name, string(year))
    @info "Exporting to $export_dir"
    raw_path = joinpath(export_dir, "$scenario_name.raw")
    md_path = joinpath(export_dir, "$(scenario_name)_metadata.json")

    # Build export files in buffers
    raw = IOBuffer()
    md = OrderedDict()
    sys = exporter.system
    psse_version = exporter.psse_version
    # These mappings are accessed in e.g. _write_bus_data via the metadata
    md["area_mapping"] = _psse_container_numbers(
        sort!(collect(PSY.get_name.(PSY.get_components(PSY.LoadZone, sys)))),
    )
    md["zone_number_mapping"] = _psse_container_numbers(
        sort!(collect(PSY.get_name.(PSY.get_components(PSY.LoadZone, sys)))),
    )
    md["record_groups"] = OrderedDict{String, Bool}()  # Keep track of which record groups we actually write to and which we skip

    # Each of these corresponds to a group of records in the PSS/E spec
    _write_case_identification_data(raw, md, sys, psse_version, "$(scenario_name)_$(year)")
    _write_bus_data(raw, md, sys, psse_version)
    _write_load_data(raw, md, sys, psse_version)
    _write_fixed_bus_shunt_data(raw, md, sys, psse_version)
    _write_generator_data(raw, md, sys, psse_version; sources_as_generators = true)
    _write_non_transformer_branch_data(raw, md, sys, psse_version)
    _write_transformer_data(raw, md, sys, psse_version)
    # TODO we'll eventually need area interchange data
    _write_skip_group(raw, md, sys, psse_version, "Area Interchange Data")
    _write_skip_group(raw, md, sys, psse_version, "Two-Terminal DC Transmission Line Data")
    _write_skip_group(raw, md, sys, psse_version,
        "Voltage Source Converter (VSC) DC Transmission Line Data")
    _write_skip_group(raw, md, sys, psse_version, "Transformer Impedance Correction Tables")
    _write_skip_group(raw, md, sys, psse_version,
        "Multi-Terminal DC Transmission Line Data")
    _write_skip_group(raw, md, sys, psse_version, "Multi-Section Line Grouping Data")
    _write_zone_data(raw, md, sys, psse_version)
    _write_skip_group(raw, md, sys, psse_version, "Interarea Transfer Data")
    _write_skip_group(raw, md, sys, psse_version, "Owner Data")
    _write_skip_group(raw, md, sys, psse_version, "FACTS Device Data")
    # TODO we'll eventually need switched shunt data
    _write_skip_group(raw, md, sys, psse_version, "Switched Shunt Data")
    _write_skip_group(raw, md, sys, psse_version, "GNE Device Data")
    _write_skip_group(raw, md, sys, psse_version, "Induction Machine Data")
    _write_q_record(raw, md, sys, psse_version)

    skipped_groups = [k for (k, v) in md["record_groups"] if !v]
    !isempty(skipped_groups) && @warn "Skipped groups: $(join(skipped_groups, ", "))"

    # Write files
    open(file -> write(file, seekstart(raw)), raw_path; truncate = true)
    open(file -> JSON.print(file, md, 4), md_path; truncate = true)

    # Write_Sienna2PSSE(exporter.system, string(scenario_name), Int64(year);
    #     export_location = string(export_location), v33 = (exporter.psse_version == :v33))
end

# TODO remove duplication between here and Write_Sienna2PSSE
"Calculate the paths of the (raw, metadata) files that would be written by a certain call to `write_export`"
function get_psse_export_paths(
    scenario_name::AbstractString,
    year::Int,
    export_location::AbstractString,
)
    base_path = joinpath(export_location, "Raw_Export", string(scenario_name), string(year))
    raw_path = joinpath(base_path, "$scenario_name.raw")
    metadata_path = joinpath(base_path, "raw_metadata_log.json")
    return (raw_path, metadata_path)
end
