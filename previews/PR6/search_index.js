var documenterSearchIndex = {"docs":
[{"location":"code_base_developer_guide/developer/#Guidelines-for-Developers","page":"Developer Guide","title":"Guidelines for Developers","text":"","category":"section"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"In order to contribute to PowerSimulationsDynamics.jl repository please read the following sections of InfrastructureSystems.jl documentation in detail:","category":"page"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"Style Guide\nContributing Guidelines","category":"page"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"Pull requests are always welcome to fix bugs or add additional modeling capabilities.","category":"page"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"All the code contributions need to include tests with a minimum coverage of 70%","category":"page"},{"location":"api/internal/#Internal","page":"Internal API Reference","title":"Internal","text":"","category":"section"},{"location":"api/internal/","page":"Internal API Reference","title":"Internal API Reference","text":"CurrentModule = PowerFlows\nDocTestSetup  = quote\n    using PowerFlows\nend","category":"page"},{"location":"api/internal/","page":"Internal API Reference","title":"Internal API Reference","text":"Modules = [PowerFlows]\nPublic = false\nPrivate = true","category":"page"},{"location":"api/internal/#PowerFlows._get_load_data-Tuple{System, Bus}","page":"Internal API Reference","title":"PowerFlows._get_load_data","text":"Obtain total load on bus b\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows._update_branch_flow!-Tuple{System}","page":"Internal API Reference","title":"PowerFlows._update_branch_flow!","text":"Updates the flow on the branches\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_func-Tuple{ACBranch, ComplexF64, ComplexF64}","page":"Internal API Reference","title":"PowerFlows.flow_func","text":"Calculates the From - To complex power flow using external data of voltages of branch of type Line\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_func-Tuple{TapTransformer, ComplexF64, ComplexF64}","page":"Internal API Reference","title":"PowerFlows.flow_func","text":"Calculates the From - To complex power flow using external data of voltages of branch of type TapTransformer\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_func-Tuple{Transformer2W, ComplexF64, ComplexF64}","page":"Internal API Reference","title":"PowerFlows.flow_func","text":"Calculates the From - To complex power flow using external data of voltages of branch of type Transformer2W\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{ACBranch}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type Line\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{DynamicBranch}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type Line\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{TapTransformer}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type TapTransformer\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{Transformer2W}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type Transformer2W\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.write_powerflow_solution!-Tuple{System, Vector{Float64}}","page":"Internal API Reference","title":"PowerFlows.write_powerflow_solution!","text":"Updates system voltages and powers with power flow results\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.write_results-Tuple{System, Vector{Float64}}","page":"Internal API Reference","title":"PowerFlows.write_results","text":"Return power flow results in dictionary of dataframes.\n\n\n\n\n\n","category":"method"},{"location":"#PowerFlows.jl","page":"Welcome Page","title":"PowerFlows.jl","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"CurrentModule = PowerFlows","category":"page"},{"location":"#Overview","page":"Welcome Page","title":"Overview","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PowerFlows.jl is a Julia package for solving Power Flows","category":"page"},{"location":"#Installation","page":"Welcome Page","title":"Installation","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"The latest stable release of PowerFlows can be installed using the Julia package manager with","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"] add PowerFlows","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"For the current development version, \"checkout\" this package with","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"] add PowerFlows#master","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PowerFlows has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP) initiative at the U.S. Department of Energy's National Renewable Energy Laboratory (NREL)","category":"page"},{"location":"api/public/#PowerFlows","page":"Public API Reference","title":"PowerFlows","text":"","category":"section"},{"location":"api/public/","page":"Public API Reference","title":"Public API Reference","text":"CurrentModule = PowerFlows\nDocTestSetup  = quote\n    using PowerFlows\nend","category":"page"},{"location":"api/public/","page":"Public API Reference","title":"Public API Reference","text":"Modules = [PowerFlows]\nPublic = true\nPrivate = false","category":"page"},{"location":"api/public/#PowerFlows.run_powerflow!-Tuple{System}","page":"Public API Reference","title":"PowerFlows.run_powerflow!","text":"Solves a the power flow into the system and writes the solution into the relevant structs. Updates generators active and reactive power setpoints and branches active and reactive power flows (calculated in the From - To direction) (see flow_val)\n\nSupports solving using Finite Differences Method (instead of using analytic Jacobian) by setting finite_diff = true. Supports passing NLsolve kwargs in the args. By default shows the solver trace.\n\nArguments available for nlsolve:\n\nmethod : See NLSolve.jl documentation for available solvers\nxtol: norm difference in x between two successive iterates under which convergence is declared. Default: 0.0.\nftol: infinite norm of residuals under which convergence is declared. Default: 1e-8.\niterations: maximum number of iterations. Default: 1_000.\nstore_trace: should a trace of the optimization algorithm's state be stored? Default: false.\nshow_trace: should a trace of the optimization algorithm's state be shown on STDOUT? Default: false.\nextended_trace: should additifonal algorithm internals be added to the state trace? Default: false.\n\nExamples\n\nrun_powerflow!(sys)\n# Passing NLsolve arguments\nrun_powerflow!(sys, method=:newton)\n# Using Finite Differences\nrun_powerflow!(sys, finite_diff=true)\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.run_powerflow-Tuple{System}","page":"Public API Reference","title":"PowerFlows.run_powerflow","text":"Similar to run_powerflow!(sys) but does not update the system struct with results. Returns the results in a dictionary of dataframes.\n\nExamples\n\nres = run_powerflow(sys)\n# Passing NLsolve arguments\nres = run_powerflow(sys, method=:newton)\n# Using Finite Differences\nres = run_powerflow(sys, finite_diff=true)\n\n\n\n\n\n","category":"method"}]
}
