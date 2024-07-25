const PSSE_EXPORT_SUPPORTED_VERSIONS = [:v33]  # TODO add :v34

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

"Peform an export from the data contained in a `PSSEExporter` to the PSS/E file format."
function write_export(
    exporter::PSSEExporter,
    scenario_name::AbstractString,
    year::Int,
    export_location::AbstractString,
)
    Write_Sienna2PSSE(exporter.system, string(scenario_name), Int64(year);
        export_location = string(export_location), v33 = (exporter.psse_version == :v33))
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
