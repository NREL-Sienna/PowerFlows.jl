const PSSE_EXPORT_SUPPORTED_VERSIONS = [:v33, :v34]

"""
Structure to perform an export from a Sienna System plus optional updates from `PowerFlowData`

# Arguments:
  - `base_system::PSY.System`: the system to be exported. Later updates may change power flow-related values but may not fundamentally alter the system
  - `psse_version::Symbol`: the version of PSS/E to target, must be one of `PSSE_EXPORT_SUPPORTED_VERSIONS`
"""
struct PSSEExporter
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
function get_paths(
    scenario_name::AbstractString,
    year::Int,
    export_location::AbstractString,
)
    base_path = joinpath(export_location, "Raw_Export", string(scenario_name), string(year))
    raw_path = joinpath(base_path, "$scenario_name.raw")
    metadata_path = joinpath(base_path, "raw_metadata_log.json")
    return (raw_path, metadata_path)
end
