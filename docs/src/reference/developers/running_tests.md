# Running Tests

We use `ReTest.jl`. The unit tests can be run in the usual way via 

```julia
julia> ] test
```

but for running the tests interactively, one can do

```julia
using TestEnv
TestEnv.activate()
include("test/load_tests.jl")
using .PowerFlowsTests
run_tests()
```

See the InfrastructureSystems.jl documentation page ["Running Tests"](https://nrel-sienna.github.io/InfrastructureSystems.jl/stable/dev_guide/tests/) for more details.

