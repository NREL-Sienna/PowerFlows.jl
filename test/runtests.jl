# See also `load_tests.jl` for running tests interactively with ReTest.jl
using PowerFlows

include("PowerFlowsTests.jl")
PowerFlowsTests.run_tests()
