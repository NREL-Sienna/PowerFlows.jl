const MAX_INIT_RESIDUAL = 10.0
const BOUNDS_TOLERANCE = 1e-6
const MAX_REACTIVE_POWER_ITERATIONS = 10
const DEFAULT_MAX_REDISTRIBUTION_ITERATIONS = 10
const LARGE_RESIDUAL = 10 # threshold for "bad initial guess": default
# norm(residual, 1)/length(residual) > 10.

const ISAPPROX_ZERO_TOLERANCE = 1e-6

const DEFAULT_NR_MAX_ITER::Int64 = 30   # default maxIterations for the NR power flow
const DEFAULT_NR_TOL::Float64 = 1e-9 # default tolerance for the NR power flow
const DEFAULT_REFINEMENT_THRESHOLD = 5e-2 # do refinement if relative error > 5%.
const DEFAULT_REFINEMENT_MAX_ITER = 10 # how many times to try iterative refinement
const DEFAULT_REFINEMENT_EPS::Float64 = 1e-6 # when to stop iterative refinement.
const NR_SINGULAR_SCALING = 1e-6 # scaling factor in fallback method for singular Jacobian
# only used for trust region.
const DEFAULT_TRUST_REGION_ETA::Float64 = 1e-4 # if actual improvement/predicted improvement
# is < eta, then reject the step Δx, shrink the trust region, and try again.
const DEFAULT_TRUST_REGION_FACTOR::Float64 = 1.0 # controls starting size of trust region
# improvement factor cutoffs for updating size of trust region.
const HALVE_TRUST_REGION = 0.1
const MAX_DOUBLE_TRUST_REGION = 0.5
const DOUBLE_TRUST_REGION = 0.9
const DEFAULT_AUTOSCALE = false # correct for scaling of the system
# typically converges in fewer iteration without autoscaling.

const PF_MAX_LOG = 10
# only used for Levenberg-Maquardt
const BIG_λ_THRESHOLD = 100.0
# 2.0 and 3.0 , or 1.5 and 5.0 for big problems; 
const DAMPING_INCR = 1.5 # no improvement => increasing damping by this factor.
const DAMPING_DECR = 5.0 # yes improvement => decrease damping by this factor.
const DEFAULT_λ_0 = 0.000001 # starting damping factor. TODO could start smaller. Not well-scaled:
# input is mix of powers (100 MW), voltages (0.8-1.2), and angles (-π/4 to π/4).
const DEFAULT_MAX_TEST_λs = 20 # give up after increasing damping factor 20 times.

const DEFAULT_Δt_k = 0.2

const AC_PF_KW = []

const PSSE_DEFAULT_EXPORT_NAME = "export"

const BUS_VOLTAGE_MAGNITUDE_CUTOFF_MIN = 0.8
const BUS_VOLTAGE_MAGNITUDE_CUTOFF_MAX = 1.2

const TIs = Union{Int32, Int64}

# voltage validation
const DEFAULT_VALIDATE_VOLTAGES = true
const MinMax = NamedTuple{(:min, :max), Tuple{Float64, Float64}}
const DEFAULT_VALIDATION_RANGE = (min = 0.5, max = 1.5)
# const MAX_INDS_TO_PRINT = 10

const OVERWRITE_NON_CONVERGED = true # overwrite non-converged time steps with NaN values

# robust homotopy method constants
const β = 10.0^-3
const INSUFFICIENT_CHANGE_IN_X = 10^(-11)
const GRAD_ZERO = 2 * eps()
# cholesky solver specific
const VTypes = SparseArrays.CHOLMOD.VRealTypes
const ITypes = SparseArrays.CHOLMOD.ITypes

# force arc names to be unique when reporting power flow results.
const FORCE_UNIQUE_NAMES = true
