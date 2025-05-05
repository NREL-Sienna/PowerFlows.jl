using Test
using Logging
using PowerSystems
using PowerSystemCaseBuilder
using PowerNetworkMatrices
using PowerFlows
using InfrastructureSystems
using LinearAlgebra
using CSV
using DataFrames
using JSON3
using Random
using DataStructures
import SparseArrays
import SparseArrays: SparseMatrixCSC, sparse, sprandn, sprand
import Random

const IS = InfrastructureSystems
const PSB = PowerSystemCaseBuilder
const PSY = PowerSystems
const PNM = PowerNetworkMatrices
const PF = PowerFlows

import Aqua
Aqua.test_unbound_args(PowerFlows)
Aqua.test_undefined_exports(PowerFlows)
Aqua.test_ambiguities(PowerFlows)
Aqua.test_stale_deps(PowerFlows)
Aqua.test_deps_compat(PowerFlows)

test_file_dir = isempty(dirname(@__FILE__)) ? "test" : dirname(@__FILE__)
const TEST_FILES_DIR = test_file_dir
const DIFF_INF_TOLERANCE = 1e-4
const DIFF_L2_TOLERANCE = 1e-3

MAIN_DIR = dirname(@__DIR__)

include("test_utils/common.jl")
include("test_utils/psse_results_compare.jl")
include("test_utils/penalty_factors_brute_force.jl")
include("test_utils/legacy_pf.jl")

const AC_SOLVERS_TO_TEST = (
    LUACPowerFlow,
    NewtonRaphsonACPowerFlow,
    TrustRegionACPowerFlow,
    LevenbergMaquardtACPowerFlow,
)

LOG_FILE = "power-flows.log"

const DISABLED_TEST_FILES = [  # Can generate with ls -1 test | grep "test_.*.jl"
# "test_dc_powerflow.jl",
# "test_iterative_methods.jl",
# "test_klu_linear_solver_cache.jl",
# "test_levenburg_marquardt.jl",
# "test_multiperiod_ac_powerflow.jl",
# "test_multiperiod_dc_powerflow.jl",
# "test_powerflow_data.jl",
# "test_psse_export.jl",
# "test_solve_powerflow.jl",
]

LOG_LEVELS = Dict(
    "Debug" => Logging.Debug,
    "Info" => Logging.Info,
    "Warn" => Logging.Warn,
    "Error" => Logging.Error,
)

"""
Copied @includetests from https://github.com/ssfrr/TestSetExtensions.jl.
Ideally, we could import and use TestSetExtensions.  Its functionality was broken by changes
in Julia v0.7.  Refer to https://github.com/ssfrr/TestSetExtensions.jl/pull/7.
"""

"""
Includes the given test files, given as a list without their ".jl" extensions.
If none are given it will scan the directory of the calling file and include all
the julia files.
"""
macro includetests(testarg...)
    if length(testarg) == 0
        tests = []
    elseif length(testarg) == 1
        tests = testarg[1]
    else
        error("@includetests takes zero or one argument")
    end

    quote
        tests = $tests
        rootfile = @__FILE__
        if length(tests) == 0
            tests = readdir(dirname(rootfile))
            tests = filter(
                f ->
                    startswith(f, "test_") && endswith(f, ".jl") && f != basename(rootfile),
                tests,
            )
        else
            tests = map(f -> string(f, ".jl"), tests)
        end
        println()
        if !isempty(DISABLED_TEST_FILES)
            @warn("Some tests are disabled $DISABLED_TEST_FILES")
        end
        for test in tests
            test in DISABLED_TEST_FILES && continue
            print(splitext(test)[1], ": ")
            include(test)
            println()
        end
    end
end

function get_logging_level_from_env(env_name::String, default)
    level = get(ENV, env_name, default)
    return IS.get_logging_level(level)
end

function run_tests()
    logging_config_filename = get(ENV, "SIIP_LOGGING_CONFIG", nothing)
    if logging_config_filename !== nothing
        config = IS.LoggingConfiguration(logging_config_filename)
    else
        config = IS.LoggingConfiguration(;
            filename = LOG_FILE,
            file_level = Logging.Info,
            console_level = Logging.Error,
        )
    end
    console_logger = ConsoleLogger(config.console_stream, config.console_level)

    IS.open_file_logger(config.filename, config.file_level) do file_logger
        levels = (Logging.Info, Logging.Warn, Logging.Error)
        multi_logger =
            IS.MultiLogger([console_logger, file_logger], IS.LogEventTracker(levels))
        global_logger(multi_logger)

        if !isempty(config.group_levels)
            IS.set_group_levels!(multi_logger, config.group_levels)
        end

        # Testing Topological components of the schema
        @time @testset "Begin PowerFlows tests" begin
            @includetests ARGS
        end

        @test length(IS.get_log_events(multi_logger.tracker, Logging.Error)) == 0
        @info IS.report_log_summary(multi_logger)
    end
end

logger = global_logger()

try
    run_tests()
finally
    # Guarantee that the global logger is reset.
    global_logger(logger)
    nothing
end
