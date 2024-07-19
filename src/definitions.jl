const RELAXED_NLSOLVE_F_TOLERANCE = :1e-4
const STRICT_NLSOLVE_F_TOLERANCE = :1e-7
const NLSOLVE_X_TOLERANCE = :1e-9
const MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE = :1e-3
const MAX_INIT_RESIDUAL = 10.0
const MAX_NLSOLVE_INTERATIONS = 10
const BOUNDS_TOLERANCE = 1e-6
const MAX_REACTIVE_POWER_ITERATIONS = 10
const DEFAULT_MAX_REDISTRIBUTION_ITERATIONS = 10

const ISAPPROX_ZERO_TOLERANCE = 1e-6
const SYSTEM_EXPORT_TOLERANCE = 1e-10  # TODO maybe this belongs in the testing codebase

const AC_PF_KW = []

# TODO remove, this is just for testing
const DATA_DIR = joinpath(dirname(dirname(dirname(pathof(PowerFlows)))), "pf_data")
