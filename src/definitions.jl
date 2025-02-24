const MAX_INIT_RESIDUAL = 10.0
const BOUNDS_TOLERANCE = 1e-6
const MAX_REACTIVE_POWER_ITERATIONS = 10
const DEFAULT_MAX_REDISTRIBUTION_ITERATIONS = 10

const ISAPPROX_ZERO_TOLERANCE = 1e-6

const DEFAULT_NR_MAX_ITER::Int64 = 30   # default maxIterations for the NR power flow
const DEFAULT_NR_TOL::Float64 = 1e-9 # default tolerance for the NR power flow

const AC_PF_KW = []

const PSSE_DEFAULT_EXPORT_NAME = "export"

const BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN = 0.8
const BUS_VOLTAGE_MAGNITUDE_CUTOFF_MAX = 1.2
