using Test
using Logging
using PowerSystems
import PowerSystems: PowerSystemTableData
const PSY = PowerSystems

import Aqua
Aqua.test_unbound_args(PowerFlows)
Aqua.test_undefined_exports(PowerFlows)
Aqua.test_ambiguities(PowerFlows)
Aqua.test_stale_deps(PowerFlows)
Aqua.test_deps_compat(PowerFlows)
