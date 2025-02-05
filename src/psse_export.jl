const PSSE_EXPORT_SUPPORTED_VERSIONS = [:v33]
const PSSE_DEFAULT = ""  # Used below in cases where we want to insert an empty field to signify the PSSE default
const PSSE_INFINITY = 9999.0
const PSSE_BUS_TYPE_MAP = Dict(
    PSY.ACBusTypes.PQ => 1,
    PSY.ACBusTypes.PV => 2,
    PSY.ACBusTypes.REF => 3,
    PSY.ACBusTypes.SLACK => 3,
    PSY.ACBusTypes.ISOLATED => 4,
)
const PSSE_BRANCH_SPECIAL_CHARACTERS = ["&", "@", "*"]

# Each of the groups in the PSS/3 v33 standard
const PSSE_GROUPS_33 = [
    "Case Identification Data",
    "Bus Data",
    "Load Data",
    "Fixed Shunt Data",
    "Generator Data",
    "Non-Transformer Branch Data",
    "Transformer Data",
    "Area Interchange Data",
    "Two-Terminal DC Transmission Line Data",
    "Voltage Source Converter (VSC) DC Transmission Line Data",
    "Transformer Impedance Correction Tables",
    "Multi-Terminal DC Transmission Line Data",
    "Multi-Section Line Grouping Data",
    "Zone Data",
    "Interarea Transfer Data",
    "Owner Data",
    "FACTS Device Data",
    "Switched Shunt Data",
    "GNE Device Data",
    "Induction Machine Data",
    "Q Record",
]

const PSSE_RAW_BUFFER_SIZEHINT = 1024
const PSSE_MD_BUFFER_SIZEHINT = 1024

"""
Structure to perform an export from a Sienna System, plus optional updates from
`PowerFlowData`, to the PSS/E format. Construct from a `System` and a PSS/E version, update
using `update_exporter` with any new data as relevant, and perform the export with
`write_export`. Writes a `<name>.raw` file and a `<name>_export_metadata.json` file with
transformations that had to be made to conform to PSS/E naming rules, which can be parsed by
PowerSystems.jl to perform a round trip with the names restored.

# Arguments:
  - `base_system::PSY.System`: the system to be exported. Later updates may change power
    flow-related values but may not fundamentally alter the system
  - `psse_version::Symbol`: the version of PSS/E to target, must be one of
    `PSSE_EXPORT_SUPPORTED_VERSIONS`
  - `write_comments::Bool` = false: whether to add the customary-but-not-in-spec-annotations
    after a slash on the first line and at group boundaries
  - `name::AbstractString = "export"`: the base name of the export
  - `step::Any = nothing`: optional step data to append to the base export name. User is
    responsible for updating the step data. If the step data is `nothing`, it is not used;
    if it is a tuple or vector, it is joined with '_' and concatted; else it is concatted
    after '_'.
  - `overwrite::Bool = false`: `true` to silently overwrite existing exports, `false` to
    throw an error if existing results are encountered
"""
mutable struct PSSEExporter <: SystemPowerFlowContainer
    system::PSY.System
    psse_version::Symbol
    export_dir::AbstractString
    name::AbstractString
    write_comments::Bool
    overwrite::Bool
    step::Any
    raw_buffer::IOBuffer  # Persist an IOBuffer to reduce allocations on repeated exports
    md_dict::OrderedDict{String}  # Persist metadata to avoid unnecessary recomputation
    md_valid::Bool  # If this is true, the metadata need not be reserialized
    md_buffer::IOBuffer  # Cache a serialized version of the metadata
    components_cache::Dict{String}  # Cache sorted lists of components to reduce allocations

    function PSSEExporter(
        base_system::PSY.System,
        psse_version::Symbol,
        export_dir::AbstractString;
        name::AbstractString = PSSE_DEFAULT_EXPORT_NAME,
        write_comments::Bool = false,
        overwrite::Bool = false,
        step::Any = nothing,
    )
        (psse_version in PSSE_EXPORT_SUPPORTED_VERSIONS) ||
            throw(
                ArgumentError(
                    "PSS/E version $psse_version is not supported, must be one of $PSSE_EXPORT_SUPPORTED_VERSIONS",
                ),
            )
        system = PSY.fast_deepcopy_system(base_system)
        mkpath(export_dir)
        new(
            system,
            psse_version,
            export_dir,
            name,
            write_comments,
            overwrite,
            step,
            IOBuffer(),
            OrderedDict{String, Any}(),
            false,
            IOBuffer(),
            Dict{String, Any}(),
        )
    end
end

supports_multi_period(::PSSEExporter) = false

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

"Force all cached information (serialized metadata, component lists, etc.) to be regenerated"
function reset_caches(exporter::PSSEExporter)
    take!(exporter.raw_buffer)
    empty!(exporter.components_cache)
    exporter.md_valid = false
    # We do not clear the md_buffer here, but !md_valid implies that its contents are not valid
end

"""
Update the `PSSEExporter` with new `data`.

# Arguments:
  - `exporter::PSSEExporter`: the exporter to update
  - `data::PSY.System`: system containing the new data. Must be fundamentally the same
  `System` as the one with which the exporter was constructed, just with different values —
  this is the user's responsibility, we do not exhaustively verify it.
"""
function update_exporter!(exporter::PSSEExporter, data::PSY.System)
    _validate_same_system(exporter.system, data) || throw(
        ArgumentError(
            "System passed to update_exporter must be the same system as the one with which the exporter was constructed, just with different values",
        ),
    )
    exporter.system = PSY.fast_deepcopy_system(data)
    reset_caches(exporter)
    return
end

get_data_array(buf::Base.GenericIOBuffer{<:Array{UInt8}}) =  # < Julia 1.11
    buf.data

(@isdefined GenericMemory) && (  # >= Julia 1.11
    get_data_array(buf::Base.GenericIOBuffer{<:GenericMemory{:not_atomic, UInt8}}) =
        Base.wrap(Array, buf.data)
)

const _FloatToBufSupportedTypes = if (@isdefined GenericMemory)
    Union{
        Base.GenericIOBuffer{<:Array{UInt8}},
        Base.GenericIOBuffer{<:GenericMemory{:not_atomic, UInt8}},
    }
else
    Base.GenericIOBuffer{<:Array{UInt8}}
end

(IOBuffer <: _FloatToBufSupportedTypes) ||
    @warn "Fast Float64 to IOBuffer implementation is out of date, will not be used"

"Temporary, very specialized proof of concept patch for https://github.com/JuliaLang/julia/issues/55835"
function better_float_to_buf(buf::_FloatToBufSupportedTypes, n::Float64)
    Base.ensureroom(buf, Base.Ryu.neededdigits(Float64))
    # get_data_array incurs an allocation on Julia >= 1.11. I think writeshortest could work
    # with the underlying Memory with minimal modification, which would be nice because
    # other than this, better_float_to_buf is completely allocation free.
    data_array = get_data_array(buf)
    new_pos = Base.Ryu.writeshortest(data_array, buf.ptr, n, false, false, true, -1,
        UInt8('e'), false, UInt8('.'), true, false)
    buf.ptr = new_pos
    buf.size = new_pos - 1
    return
end

fastprint(io::IO, val) = print(io, val)
fastprint(io::_FloatToBufSupportedTypes, val::Float64) = better_float_to_buf(io, val)

function fastprintdelim(io, val, delim = ", ")
    fastprint(io, val)
    fastprint(io, delim)
end

function fastprintln(io, val, ln = "\n")
    fastprint(io, val)
    fastprint(io, ln)
end

# call fastprintdelim multiple times on val, unrolled
macro fastprintdelim_multi(io, val, newline::Bool, mult::Int)
    local exprs = [Expr(:call, :fastprintdelim, esc(io), esc(val)) for _ in 1:(mult - 1)]
    local lastCall = newline ? :fastprintln : :fastprintdelim
    local lastExpr = Expr(:call, lastCall, esc(io), esc(val))
    return Expr(:block, exprs..., lastExpr)
end
# call fastprintdelim on each item in vals, unrolled
macro fastprintdelim_unroll(io, newline::Bool, vals...)
    local allButLast = vals[begin:(end - 1)]
    local exprs = [Expr(:call, :fastprintdelim, esc(io), esc(item)) for item in allButLast]
    local lastCall = newline ? :fastprintln : :fastprintdelim
    local lastExpr = Expr(:call, lastCall, esc(io), esc(vals[end]))
    return Expr(:block, exprs..., lastExpr)
end

function fastprintdelim_psse_default_ownership(io)
    @fastprintdelim_multi(io, PSSE_DEFAULT, false, 8)
end

function fastprintln_psse_default_ownership(io)
    @fastprintdelim_multi(io, PSSE_DEFAULT, true, 8)
end

function end_group_33(io::IO, md::AbstractDict, exporter::PSSEExporter, group_name, written)
    next_group = PSSE_GROUPS_33[only(findall(==(group_name), PSSE_GROUPS_33)) + 1]
    end_msg = "0"
    if exporter.write_comments
        end_msg *= " / End of $group_name"
        (next_group == "Q Record") || (end_msg *= ", Begin $next_group")
    end
    println(io, end_msg)
    exporter.md_valid || (md["record_groups"][group_name] = written)
end

# Parses "1" and "1.0" as 1, returns nothing on "1.5" and "a"
function _permissive_parse_int(x)
    n = tryparse(Float64, x)
    isnothing(n) && return nothing
    (round(n) == n) || return nothing
    return Int64(n)
end

"""
If `val` is empty, returns `T()`; if not, asserts that `val isa T` and returns `val`. Has nice type checker semantics.

# Examples
```julia
convert_empty(Vector{String}, [])  # -> String[]
convert_empty(Vector{String}, ["a"])  # -> ["a"]
convert_empty(Vector{String}, [2])  # -> TypeError: in typeassert, expected Vector{String}, got a value of type Vector{Int64}
Base.return_types(Base.Fix1(convert_empty, Vector{String}))  # -> [Vector{String}]
```
"""
convert_empty(::Type{T}, val) where {T} = isempty(val) ? T() : val::T
convert_empty_stringvec = Base.Fix1(convert_empty, Vector{String})

# PERF could be improved by appending to the buffer rather than doing string interpolation, seems unnecessary
_psse_quote_string(s::String) = "'$s'"

branch_to_bus_numbers(branch) =
    PSY.get_number.((PSY.get_from_bus(branch), PSY.get_to_bus(branch)))::Tuple{Int, Int}

"Throw a `NotImplementedError` if the `psse_version` is not `:v33`"
check_33(exporter::PSSEExporter) = check_33(exporter.psse_version)
check_33(psse_version::Symbol) =
    (psse_version == :v33) ||
    throw(IS.NotImplementedError("Only implemented for psse_version $(:v33)"))

"Validate that the Sienna area/zone names parse as PSS/E-compatible area/zone numbers, output a mapping"
function _map_psse_container_names(container_names::Vector{String})
    (length(container_names) <= 9999) || throw(ArgumentError("Too many container_names"))
    mapping = OrderedDict{String, Int}()
    used = Set{Int}()
    for name in container_names
        try
            # Ideally it's a number and we just parse it
            parsed = _permissive_parse_int(name)
            if !isnothing(parsed) && (parsed in 1:9999) && !(parsed in used)
                mapping[name] = parsed
                push!(used, parsed)
                continue
            end
        catch e
            (e isa Union{InexactError, ArgumentError}) || rethrow(e)
        end

        # If parsing fails or the number doesn't work, assign it the lowest unused number
        i = 1
        while i in used  # PERF inefficient but unlikely to matter
            i += 1
        end
        @assert i <= 9999
        mapping[name] = i
        push!(used, i)
    end
    return mapping
end

"WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Case Identification Data"
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Case Identification Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict

    check_33(exporter)
    now = Dates.now()
    md_string = "PSS/E 33.3 RAW via PowerFlows.jl, $now"

    # Record 1
    IC = 0
    SBASE = PSY.get_base_power(exporter.system)
    REV = 33
    XFRRAT = 0
    NXFRAT = 1
    BASFRQ = PSY.get_frequency(exporter.system)
    exporter.write_comments && (BASFRQ = "$BASFRQ    / $md_string")

    @fastprintdelim_unroll(io, true, IC, SBASE, REV, XFRRAT, NXFRAT, BASFRQ)

    # Record 2
    case_name = md["case_name"]
    (length(case_name) <= 60) ||
        throw(ArgumentError("case_name may be up to 60 characters"))
    println(io, case_name)

    # Record 3
    line3 = md_string
    @assert length(line3) <= 60
    println(io, line3)
    exporter.md_valid || (md["record_groups"]["Case Identification Data"] = true)
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
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Bus Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)

    buses = get!(exporter.components_cache, "buses") do
        sort!(collect(PSY.get_components(PSY.Bus, exporter.system)); by = PSY.get_number)
    end
    old_bus_numbers = PSY.get_number.(buses)

    if !exporter.md_valid
        md["bus_number_mapping"] = _psse_bus_numbers(old_bus_numbers)
        md["bus_name_mapping"] =
            _psse_bus_names(
                convert_empty_stringvec(PSY.get_name.(buses)),
                old_bus_numbers,
                md["bus_number_mapping"],
            )
    end
    bus_number_mapping = md["bus_number_mapping"]
    bus_name_mapping = md["bus_name_mapping"]

    for bus in buses
        I = bus_number_mapping[PSY.get_number(bus)]
        NAME = _psse_quote_string(bus_name_mapping[PSY.get_name(bus)])
        BASKV = PSY.get_base_voltage(bus)
        IDE = PSSE_BUS_TYPE_MAP[PSY.get_bustype(bus)]
        AREA = if isnothing(PSY.get_area(bus))
            PSSE_DEFAULT
        else
            md["area_mapping"][PSY.get_name(PSY.get_area(bus))]
        end
        ZONE = if isnothing(PSY.get_load_zone(bus))
            PSSE_DEFAULT
        else
            md["zone_mapping"][PSY.get_name(PSY.get_load_zone(bus))]
        end
        OWNER = PSSE_DEFAULT
        VM = PSY.get_magnitude(bus)
        VA = rad2deg(PSY.get_angle(bus))
        NVHI = PSY.get_voltage_limits(bus).max
        NVLO = PSY.get_voltage_limits(bus).min
        EVHI = PSSE_DEFAULT
        EVLO = PSSE_DEFAULT

        @fastprintdelim_unroll(io, true, I, NAME, BASKV, IDE, AREA,
            ZONE, OWNER, VM, VA,
            NVHI, NVLO, EVHI, EVLO)
    end
    end_group_33(io, md, exporter, "Bus Data", true)
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

"Take the output of `create_component_ids` and make it more suitable for JSON serialization"
serialize_component_ids(id_mapping::Dict{Tuple{Int64, String}, String}) =
    Dict("$(s_bus_n)_$(s_name)" => p_name for ((s_bus_n, s_name), p_name) in id_mapping)
serialize_component_ids(id_mapping::Dict{Tuple{Tuple{Int64, Int64}, String}, String}) =
    Dict(
        "$(s_bus_1)-$(s_bus_2)_$(s_name)" => p_name for
        (((s_bus_1, s_bus_2), s_name), p_name) in id_mapping
    )

# Fetch PL, QL, IP, IQ, YP, YQ
_psse_get_load_data(exporter::PSSEExporter, load::PSY.StandardLoad) =
    with_units_base(exporter.system, PSY.UnitSystem.NATURAL_UNITS) do
        PSY.get_constant_active_power(load),
        PSY.get_constant_reactive_power(load),
        PSY.get_current_active_power(load),
        PSY.get_current_reactive_power(load),
        PSY.get_impedance_active_power(load),
        PSY.get_impedance_reactive_power(load)
    end

# Fallback if not all the data is available
# This mapping corresponds to `function make_power_load` in the parser
_psse_get_load_data(exporter::PSSEExporter, load::PSY.StaticLoad) =
    with_units_base(exporter.system, PSY.UnitSystem.NATURAL_UNITS) do
        PSY.get_active_power(load),
        PSY.get_reactive_power(load),
        PSSE_DEFAULT,
        PSSE_DEFAULT,
        PSSE_DEFAULT,
        PSSE_DEFAULT
    end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Load Data
"""
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Load Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)

    loads = get!(exporter.components_cache, "loads") do
        sort!(collect(PSY.get_components(PSY.StaticLoad, exporter.system)); by = PSY.get_name)
    end
    load_name_mapping = get!(exporter.components_cache, "load_name_mapping") do
        create_component_ids(
            convert_empty_stringvec(PSY.get_name.(loads)),
            PSY.get_number.(PSY.get_bus.(loads));
            singles_to_1 = true,
        )
    end
    for load in loads
        sienna_bus_number = PSY.get_number(PSY.get_bus(load))
        I = md["bus_number_mapping"][sienna_bus_number]
        ID = _psse_quote_string(load_name_mapping[(sienna_bus_number, PSY.get_name(load))])
        STATUS = PSY.get_available(load) ? 1 : 0
        AREA = PSSE_DEFAULT  # defaults to bus's area
        ZONE = PSSE_DEFAULT  # defaults to zone's area
        PL, QL, IP, IQ, YP, YQ = _psse_get_load_data(exporter, load)
        OWNER = PSSE_DEFAULT  # defaults to bus's owner
        SCALE = PSSE_DEFAULT
        INTRPT = PSSE_DEFAULT

        @fastprintdelim_unroll(io, true, I, ID, STATUS, AREA, ZONE,
            PL, QL, IP, IQ, YP, YQ, OWNER,
            SCALE, INTRPT)
    end
    end_group_33(io, md, exporter, "Load Data", true)
    exporter.md_valid ||
        (md["load_name_mapping"] = serialize_component_ids(load_name_mapping))
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Fixed Bus Shunt Data
"""
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Fixed Shunt Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)

    shunts = get!(exporter.components_cache, "shunts") do
        sort!(
            collect(PSY.get_components(PSY.FixedAdmittance, exporter.system));
            by = PSY.get_name,
        )
    end
    shunt_name_mapping = get!(exporter.components_cache, "shunt_name_mapping") do
        create_component_ids(
            convert_empty_stringvec(PSY.get_name.(shunts)),
            PSY.get_number.(PSY.get_bus.(shunts));
            singles_to_1 = true,
        )
    end
    for shunt in shunts
        sienna_bus_number = PSY.get_number(PSY.get_bus(shunt))
        I = md["bus_number_mapping"][sienna_bus_number]
        ID =
            _psse_quote_string(shunt_name_mapping[(sienna_bus_number, PSY.get_name(shunt))])
        STATUS = PSY.get_available(shunt) ? 1 : 0
        GL = real(PSY.get_Y(shunt)) * PSY.get_base_power(exporter.system)
        BL = imag(PSY.get_Y(shunt)) * PSY.get_base_power(exporter.system)

        @fastprintdelim_unroll(io, true, I, ID, STATUS, GL, BL)
    end
    end_group_33(io, md, exporter, "Fixed Shunt Data", true)
    exporter.md_valid ||
        (md["shunt_name_mapping"] = serialize_component_ids(shunt_name_mapping))
end

function _warn_finite_default(val; field_name, component_name)
    isfinite(val) && return val
    if val == Inf
        newval = PSSE_INFINITY
    elseif val == -Inf
        newval = -PSSE_INFINITY
    elseif isnan(val)
        newval = PSSE_DEFAULT
    else
        error("Should be unreachable")
    end
    @warn "Detected non-finite value $field_name = $val for $component_name, using '$newval'"
    return newval
end

"""
If the export_settings flag `sources_as_generators` is set, export `PSY.Source` instances as
PSS/E generators in addition to `PSY.Generator`s. Same for `storages_as_generators` and
`PSY.Storage`.

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Generator Data
"""
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Generator Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)

    generators = get!(exporter.components_cache, "generators") do
        temp_gens::Vector{PSY.StaticInjection} = sort!(
            # TODO consider storage here
            collect(PSY.get_components(PSY.Generator, exporter.system));
            by = PSY.get_name,
        )
        get(md["export_settings"], "sources_as_generators", false) && append!(
            temp_gens,
            sort!(
                collect(PSY.get_components(PSY.Source, exporter.system));
                by = PSY.get_name,
            ),
        )
        get(md["export_settings"], "storages_as_generators", false) && append!(
            temp_gens,
            sort!(
                collect(PSY.get_components(PSY.Storage, exporter.system));
                by = PSY.get_name,
            ),
        )
        return temp_gens
    end
    generator_name_mapping = get!(exporter.components_cache, "generator_name_mapping") do
        create_component_ids(
            convert_empty_stringvec(PSY.get_name.(generators)),
            PSY.get_number.(PSY.get_bus.(generators));
            singles_to_1 = false,
        )
    end
    for generator in generators
        sienna_bus_number = PSY.get_number(PSY.get_bus(generator))
        I = md["bus_number_mapping"][sienna_bus_number]
        ID =
            _psse_quote_string(
                generator_name_mapping[(sienna_bus_number, PSY.get_name(generator))],
            )
        PG, QG = with_units_base(exporter.system, PSY.UnitSystem.NATURAL_UNITS) do
            PSY.get_active_power(generator), PSY.get_reactive_power(generator)
        end
        reactive_power_limits = with_units_base(
            () -> get_reactive_power_limits_for_power_flow(generator),
            exporter.system,
            PSY.UnitSystem.NATURAL_UNITS,
        )
        QT = reactive_power_limits.max
        QT = _warn_finite_default(
            QT;
            field_name = "QT",
            component_name = PSY.get_name(generator),
        )
        QB = reactive_power_limits.min
        QB = _warn_finite_default(
            QB;
            field_name = "QB",
            component_name = PSY.get_name(generator),
        )
        VS = PSY.get_magnitude(PSY.get_bus(generator))
        IREG = get(PSY.get_ext(generator), "IREG", PSSE_DEFAULT)
        MBASE = PSY.get_base_power(generator)
        ZR, ZX = PSSE_DEFAULT, PSSE_DEFAULT
        RT, XT = PSSE_DEFAULT, PSSE_DEFAULT
        GTAP = PSSE_DEFAULT
        STAT = PSY.get_available(generator) ? 1 : 0
        RMPCT = PSSE_DEFAULT
        active_power_limits =
            with_units_base(
                () -> get_active_power_limits_for_power_flow(generator),
                exporter.system,
                PSY.UnitSystem.NATURAL_UNITS,
            )
        PT = active_power_limits.max
        PT = _warn_finite_default(
            PT;
            field_name = "PT",
            component_name = PSY.get_name(generator),
        )
        PB = active_power_limits.min
        PB = _warn_finite_default(
            PB;
            field_name = "PB",
            component_name = PSY.get_name(generator),
        )
        WMOD = get(PSY.get_ext(generator), "WMOD", PSSE_DEFAULT)
        WPF = get(PSY.get_ext(generator), "WPF", PSSE_DEFAULT)

        @fastprintdelim_unroll(io, false, I, ID, PG, QG, QT, QB,
            VS, IREG, MBASE, ZR, ZX,
            RT, XT, GTAP, STAT, RMPCT,
            PT, PB)
        fastprintdelim_psse_default_ownership(io)
        @fastprintdelim_unroll(io, true, WMOD, WPF)
    end
    end_group_33(io, md, exporter, "Generator Data", true)
    exporter.md_valid ||
        (md["generator_name_mapping"] = serialize_component_ids(generator_name_mapping))
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Non-Transformer Branch Data
"""
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Non-Transformer Branch Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)

    branches_with_numbers = get!(exporter.components_cache, "branches") do
        branches = sort!(
            collect(
                PSY.get_components(Union{PSY.Line, PSY.MonitoredLine}, exporter.system),
            );
            by = branch_to_bus_numbers,
        )
        [(branch, branch_to_bus_numbers(branch)) for branch in branches]
    end
    branch_name_mapping = get!(exporter.components_cache, "branch_name_mapping") do
        create_component_ids(
            convert_empty_stringvec(PSY.get_name.(first.(branches_with_numbers))),
            last.(branches_with_numbers);
            singles_to_1 = false,
        )
    end

    for (branch, (from_n, to_n)) in branches_with_numbers
        I = md["bus_number_mapping"][from_n]
        J = md["bus_number_mapping"][to_n]
        CKT = branch_name_mapping[((from_n, to_n), PSY.get_name(branch))]
        if first(CKT) in PSSE_BRANCH_SPECIAL_CHARACTERS
            if first(CKT) == '&'
                @error "Exporting branch $(PSY.get_name(branch)) with disallowed name '$CKT'"
            elseif first(CKT) == '@'
                @warn "Exporting branch $(PSY.get_name(branch)) with name '$CKT', it will be treated as a breaker"
            elseif first(CKT) == '*'
                @warn "Exporting branch $(PSY.get_name(branch)) with name '$CKT', it will be treated as a switch"
            else
                error("Should be unreachable")
            end
        end
        # '&', 
        CKT = _psse_quote_string(CKT)
        R = PSY.get_r(branch)
        X = PSY.get_x(branch)
        B = 0.0  # NOTE PowerSystems only represents BI, BJ
        RATEA =
            RATEB =
                RATEC =
                    with_units_base(
                        () -> PSY.get_rating(branch),
                        exporter.system,
                        PSY.UnitSystem.NATURAL_UNITS,
                    )
        GI, BI = 0.0, PSY.get_b(branch).from
        GJ, BJ = 0.0, PSY.get_b(branch).to
        ST = PSY.get_available(branch) ? 1 : 0
        MET = PSSE_DEFAULT
        LEN = PSSE_DEFAULT

        @fastprintdelim_unroll(io, false, I, J, CKT, R, X, B,
            RATEA, RATEB, RATEC, GI, BI,
            GJ, BJ, ST, MET, LEN)
        fastprintln_psse_default_ownership(io)
    end
    end_group_33(io, md, exporter, "Non-Transformer Branch Data", true)
    exporter.md_valid ||
        (md["branch_name_mapping"] = serialize_component_ids(branch_name_mapping))
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
        # If both bus numbers are large, that new_name might be too long
        if !_is_valid_psse_name(new_name)
            n = 0
            while !_is_valid_psse_name(new_name) || (new_name in used_names)
                new_name = "B$(bus_number_mapping[orig_from])-N$n"
                n += 1
            end
        end
        @assert _is_valid_psse_name(new_name) new_name
        mapping[original_name] = new_name
        push!(used_names, new_name)
    end
    return mapping
end

"""
Currently only supports two-winding transformers

WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Transformer Data
"""
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Transformer Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)
    transformer_types =
        Union{PSY.Transformer2W, PSY.TapTransformer, PSY.PhaseShiftingTransformer}
    transformers_with_numbers = get!(exporter.components_cache, "transformers") do
        transformers = sort!(
            collect(PSY.get_components(transformer_types, exporter.system));
            by = branch_to_bus_numbers,
        )
        [(transformer, branch_to_bus_numbers(transformer)) for transformer in transformers]
    end
    transformer_ckt_mapping = get!(exporter.components_cache, "transformer_ckt_mapping") do
        create_component_ids(
            convert_empty_stringvec(PSY.get_name.(first.(transformers_with_numbers))),
            last.(transformers_with_numbers);
            singles_to_1 = false,
        )
    end
    if !exporter.md_valid
        md["transformer_name_mapping"] = _psse_transformer_names(
            convert_empty_stringvec(PSY.get_name.(first.(transformers_with_numbers))),
            last.(transformers_with_numbers),
            md["bus_number_mapping"],
            transformer_ckt_mapping,
        )
    end

    bus_number_mapping = md["bus_number_mapping"]
    transformer_name_mapping = md["transformer_name_mapping"]

    for (transformer, (from_n, to_n)) in transformers_with_numbers
        from_n, to_n = branch_to_bus_numbers(transformer)
        I = bus_number_mapping[from_n]
        J = bus_number_mapping[to_n]
        K = 0  # no third winding
        CKT = transformer_ckt_mapping[((from_n, to_n), PSY.get_name(transformer))]
        @assert !(first(CKT) in ['&', '@', '*'])  # Characters with a special meaning in this context
        CKT = _psse_quote_string(CKT)
        CW = 1  # NOTE on parsing we do the transformation to this unit system
        CZ = 1  # NOTE on parsing we do the transformation to this unit system
        CM = 1  # NOTE on parsing we do the transformation to this unit system
        MAG1 = PSSE_DEFAULT
        MAG2 = -PSY.get_primary_shunt(transformer)
        (MAG2 > 0) &&
            @warn "Detected positive MAG2 for transformer $(PSY.get_name(transformer)) due to negative `get_primary_shunt`; exporting anyway"
        NMETR = PSSE_DEFAULT
        NAME = _psse_quote_string(transformer_name_mapping[PSY.get_name(transformer)])
        STAT = PSY.get_available(transformer) ? 1 : 0
        VECGRP = PSSE_DEFAULT

        R1_2 = PSY.get_r(transformer)
        X1_2 = PSY.get_x(transformer)
        SBASE1_2 = PSY.get_base_power(transformer)

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
                    with_units_base(
                        () -> PSY.get_rating(transformer),
                        exporter.system,
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

        @fastprintdelim_unroll(io, false, I, J, K, CKT, CW, CZ, CM,
            MAG1, MAG2, NMETR, NAME, STAT)
        fastprintdelim_psse_default_ownership(io)
        fastprintln(io, VECGRP)

        @fastprintdelim_unroll(io, true, R1_2, X1_2, SBASE1_2)

        @fastprintdelim_unroll(io, true, WINDV1, NOMV1, ANG1, RATA1,
            RATB1, RATC1, COD1, CONT1, RMA1, RMI1,
            VMA1, VMI1, NTP1, TAB1, CR1, CX1, CNXA1)

        @fastprintdelim_unroll(io, true, WINDV2, NOMV2)
    end
    end_group_33(io, md, exporter, "Transformer Data", true)
    if !exporter.md_valid
        md["transformer_ckt_mapping"] = serialize_component_ids(transformer_ckt_mapping)
    end
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Zone Data
"""
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Zone Data")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)
    zone_mapping = md["zone_mapping"]
    zones = get!(exporter.components_cache, "zones") do
        sort!(
            collect(PSY.get_components(PSY.LoadZone, exporter.system));
            by = x -> zone_mapping[PSY.get_name(x)],
        )
    end
    for zone in zones
        name = PSY.get_name(zone)
        I = zone_mapping[name]
        @assert _is_valid_psse_name(name) name
        ZONAME = _psse_quote_string(name)

        @fastprintdelim_unroll(io, true, I, ZONAME)
    end
    end_group_33(io, md, exporter, "Zone Data", true)
end

"""
WRITTEN TO SPEC: PSS/E 33.3 POM 5.2.1 Q Record
"""
function write_to_buffers!(
    exporter::PSSEExporter,
    ::Val{Symbol("Q Record")},
)
    io = exporter.raw_buffer
    md = exporter.md_dict
    check_33(exporter)
    println(io, "Q")  # End of file
    exporter.md_valid || (md["record_groups"]["Q Record"] = true)
end

function _write_skip_group(
    io::IO,
    md::AbstractDict,
    exporter::PSSEExporter,
    this_section_name::String,
)
    check_33(exporter)
    end_group_33(io, md, exporter, this_section_name, false)
    exporter.md_valid || (md["record_groups"][this_section_name] = false)
end

# If a writer for a given group is not defined, write that we are skipping it
function write_to_buffers!(exporter::PSSEExporter, ::Val{T}) where {T}
    io = exporter.raw_buffer
    md = exporter.md_dict
    group_name = string(T)
    @debug "Export for group $group_name not implemented, skipping it"
    _write_skip_group(io, md, exporter, group_name)
end

_step_to_string(::Nothing) = ""
_step_to_string(iterable_step::Union{Tuple, AbstractArray}) = "_" * join(iterable_step, "_")
_step_to_string(scalar_step::Any) = "_$scalar_step"

"Peform an export from the data contained in a `PSSEExporter` to the PSS/E file format."
function write_export(
    exporter::PSSEExporter,
    name::AbstractString;
    overwrite = false,
)
    original_name = name
    name = name * _step_to_string(exporter.step)
    # Construct paths
    export_subdir = joinpath(exporter.export_dir, name)
    dir_exists = isdir(export_subdir)
    (dir_exists && !overwrite) && throw(
        ArgumentError(
            "Target export directory $(abspath(export_subdir)) already exists; specify `overwrite = true` if it should be overwritten",
        ),
    )
    dir_exists || mkdir(export_subdir)
    @info "Exporting $name to $export_subdir"
    raw_path, md_path = get_psse_export_paths(export_subdir)

    # Build export files in buffers
    md = exporter.md_dict
    if !exporter.md_valid
        md["case_name"] = name

        md["export_settings"] = OrderedDict{String, Any}()
        export_settings = md["export_settings"]
        export_settings["psse_version"] = string(exporter.psse_version)
        export_settings["export_dir"] = exporter.export_dir
        export_settings["original_name"] = original_name
        export_settings["write_comments"] = exporter.write_comments
        export_settings["overwrite"] = exporter.overwrite
        export_settings["step"] = _step_to_string(exporter.step)
        export_settings["sources_as_generators"] = true
        export_settings["storages_as_generators"] = true

        md["record_groups"] = OrderedDict{String, Bool}()  # Keep track of which record groups we actually write to and which we skip

        # These mappings are accessed in e.g. _write_bus_data via the metadata
        md["area_mapping"] = _map_psse_container_names(
            sort!(
                collect(
                    convert_empty_stringvec(
                        PSY.get_name.(PSY.get_components(PSY.Area, exporter.system)),
                    ),
                )),
        )
        md["zone_mapping"] = _map_psse_container_names(
            sort!(
                collect(
                    convert_empty_stringvec(
                        PSY.get_name.(PSY.get_components(PSY.LoadZone, exporter.system)),
                    ),
                )),
        )
    end

    with_units_base(exporter.system, PSY.UnitSystem.SYSTEM_BASE) do
        # Each of these corresponds to a group of records in the PSS/E spec
        for group_name in PSSE_GROUPS_33
            @debug "Writing export for group $group_name"
            write_to_buffers!(exporter, Val{Symbol(group_name)}())
        end
    end

    skipped_groups = [k for (k, v) in md["record_groups"] if !v]
    !isempty(skipped_groups) && @warn "Skipped groups: $(join(skipped_groups, ", "))"

    if !exporter.md_valid
        @debug "Serializing PSS/E export metadata to in-memory buffer"
        take!(exporter.md_buffer)
        JSON3.pretty(exporter.md_buffer, md)
    end
    exporter.md_valid = true

    # Write files
    open(file -> write(file, take!(exporter.raw_buffer)), raw_path; truncate = true)
    open(file -> write(file, seekstart(exporter.md_buffer)), md_path; truncate = true)
    return
end

write_export(exporter::PSSEExporter; kwargs...) =
    write_export(
        exporter,
        exporter.name;
        merge(Dict(:overwrite => exporter.overwrite), kwargs)...,
    )

"Calculate the paths of the (raw, metadata) files that would be written by a certain call to `write_export`"
function get_psse_export_paths(
    export_subdir::AbstractString,
)
    name = last(splitdir(export_subdir))
    raw_path = joinpath(export_subdir, "$name.raw")
    metadata_path = joinpath(export_subdir, "$(name)$(PSY.PSSE_EXPORT_METADATA_EXTENSION)")
    return (raw_path, metadata_path)
end

# COMMON INTERFACE
make_power_flow_container(pfem::PSSEExportPowerFlow, sys::PSY.System; kwargs...) =
    PSSEExporter(
        sys,
        pfem.psse_version,
        pfem.export_dir;
        name = pfem.name,
        write_comments = pfem.write_comments,
        overwrite = pfem.overwrite,
        step = (0, 0),
    )

solve_powerflow!(exporter::PSSEExporter) = write_export(exporter)
