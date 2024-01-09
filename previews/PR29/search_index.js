var documenterSearchIndex = {"docs":
[{"location":"code_base_developer_guide/developer/#Guidelines-for-Developers","page":"Developer Guide","title":"Guidelines for Developers","text":"","category":"section"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"In order to contribute to PowerSimulationsDynamics.jl repository please read the following sections of InfrastructureSystems.jl documentation in detail:","category":"page"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"Style Guide\nContributing Guidelines","category":"page"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"Pull requests are always welcome to fix bugs or add additional modeling capabilities.","category":"page"},{"location":"code_base_developer_guide/developer/","page":"Developer Guide","title":"Developer Guide","text":"All the code contributions need to include tests with a minimum coverage of 70%","category":"page"},{"location":"api/internal/#Internal","page":"Internal API Reference","title":"Internal","text":"","category":"section"},{"location":"api/internal/","page":"Internal API Reference","title":"Internal API Reference","text":"CurrentModule = PowerFlows\nDocTestSetup  = quote\n    using PowerFlows\nend","category":"page"},{"location":"api/internal/","page":"Internal API Reference","title":"Internal API Reference","text":"Modules = [PowerFlows]\nPublic = false\nPrivate = true","category":"page"},{"location":"api/internal/#PowerFlows._get_load_data-Tuple{System, Bus}","page":"Internal API Reference","title":"PowerFlows._get_load_data","text":"Obtain total load on bus b\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows._update_branch_flow!-Tuple{System}","page":"Internal API Reference","title":"PowerFlows._update_branch_flow!","text":"Updates the flow on the branches\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_func-Tuple{ACBranch, ComplexF64, ComplexF64}","page":"Internal API Reference","title":"PowerFlows.flow_func","text":"Calculates the From - To complex power flow using external data of voltages of branch of type Line\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_func-Tuple{TapTransformer, ComplexF64, ComplexF64}","page":"Internal API Reference","title":"PowerFlows.flow_func","text":"Calculates the From - To complex power flow using external data of voltages of branch of type TapTransformer\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_func-Tuple{Transformer2W, ComplexF64, ComplexF64}","page":"Internal API Reference","title":"PowerFlows.flow_func","text":"Calculates the From - To complex power flow using external data of voltages of branch of type Transformer2W\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{ACBranch}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type Line\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{DynamicBranch}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type Line\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{TapTransformer}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type TapTransformer\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.flow_val-Tuple{Transformer2W}","page":"Internal API Reference","title":"PowerFlows.flow_val","text":"Calculates the From - To complex power flow (Flow injected at the bus) of branch of type Transformer2W\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.solve_powerflow!-Tuple{PowerFlowData{PowerNetworkMatrices.PTDF{Tuple{Vector{Int64}, Vector{String}}, Tuple{Dict{Int64, Int64}, Dict{String, Int64}}, Matrix{Float64}}, PowerNetworkMatrices.ABA_Matrix{Tuple{Vector{Int64}, Vector{Int64}}, Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}}, KLU.KLUFactorization{Float64, Int64}}}}","page":"Internal API Reference","title":"PowerFlows.solve_powerflow!","text":"Evaluates the power flows on each system's branch and updates the PowerFlowData structure.\n\nArguments:\n\ndata::PTDFPowerFlowData:       PTDFPowerFlowData structure containig all the information related to the system's power flow.\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.solve_powerflow!-Tuple{PowerFlowData}","page":"Internal API Reference","title":"PowerFlows.solve_powerflow!","text":"Evaluates the power flows on each system's branch and updates the PowerFlowData structure.\n\nArguments:\n\ndata::vPTDFPowerFlowData:       vPTDFPowerFlowData structure containig all the information related to the system's power flow.\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerFlows.write_powerflow_solution!-Tuple{System, Vector{Float64}}","page":"Internal API Reference","title":"PowerFlows.write_powerflow_solution!","text":"Updates system voltages and powers with power flow results\n\n\n\n\n\n","category":"method"},{"location":"modeler_guide/power_flow/#Power-Flow","page":"Power Flow","title":"Power Flow","text":"","category":"section"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"CurrentModule = PowerFlows","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"PowerFlows.jl provides the capability to run a power flow using NLSolve, in the current stage of development it can't force reactive power constraints. This power flow routine does not check for reactive power limits or other limiting mechanisms in the grid, and can therefore be used to check for solver convergence - making no guarantees of the solution feasibility.","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"The power flow solver uses NLsolve.jl under the hood and takes any keyword argument accepted by NLsolve. The solver uses the current operating point in the buses to provide the initial guess.","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"Limitations: The PowerFlow solver doesn't support systems with HVDC lines or Phase Shifting transformers yet. The power flow solver can't handle systems with islands.","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"using PowerFlows\nusing PowerSystems\nusing PowerSystemCaseBuilder\n\nsystem_data = build_system(PSITestSystems, \"c_sys14\")","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"PowerFlows.jl has two modes of using the power flow solver.","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"Solving the power flow for the current operating point in the system. Takes the data in the buses, the active_power and reactive_power fields in the static injection devices. Returns a dictionary with results in a DataFrame that can be exported or manipulated as needed.\nSolves the power flow and updated the devices in the system to the operating condition. This model will update the values of magnitudes and angles in the system's buses. It also updates the active and reactive power flows in the branches and devices connected to PV buses. It also updates the active and reactive power of the injection devices connected to the Slack bus, and updates only the reactive power of the injection devices connected to PV buses. If multiple devices are connected to the same bus, the power is divided proportional to the base power. This utility is useful to initialize systems before serializing or checking the addition of new devices is still AC feasible.","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"Solving the power flow with mode 1:","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"results = run_powerflow(system_data)\nresults[\"bus_results\"]","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"Solving the power flow with mode 2:","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"Before running the power flow command these are the values of the voltages:","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"for b in get_components(Bus, system_data)\n    println(\"$(get_name(b)) - Magnitude $(get_magnitude(b)) - Angle (rad) $(get_angle(b))\")\nend","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"run_powerflow! return true or false to signal the successful result of the power flow. This enables the integration of a power flow into functions and use the return as check. For instance, initializing dynamic simulations. Also, because run_powerflow! uses NLsolve.jl all the parameters used for NLsolve are also available for run_powerflow!","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"run_powerflow!(system_data; method = :newton)","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"After running the power flow command this are the values of the voltages:","category":"page"},{"location":"modeler_guide/power_flow/","page":"Power Flow","title":"Power Flow","text":"for b in get_components(Bus, system_data)\n    println(\"$(get_name(b)) - Magnitude $(get_magnitude(b)) - Angle (rad) $(get_angle(b))\")\nend","category":"page"},{"location":"#PowerFlows.jl","page":"Welcome Page","title":"PowerFlows.jl","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"CurrentModule = PowerFlows","category":"page"},{"location":"#Overview","page":"Welcome Page","title":"Overview","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PowerFlows.jl is a Julia package for solving Power Flows","category":"page"},{"location":"#Installation","page":"Welcome Page","title":"Installation","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"The latest stable release of PowerFlows can be installed using the Julia package manager with","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"] add PowerFlows","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"For the current development version, \"checkout\" this package with","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"] add PowerFlows#main","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PowerFlows has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP) initiative at the U.S. Department of Energy's National Renewable Energy Laboratory (NREL)","category":"page"},{"location":"api/public/#PowerFlows","page":"Public API Reference","title":"PowerFlows","text":"","category":"section"},{"location":"api/public/","page":"Public API Reference","title":"Public API Reference","text":"CurrentModule = PowerFlows\nDocTestSetup  = quote\n    using PowerFlows\nend","category":"page"},{"location":"api/public/","page":"Public API Reference","title":"Public API Reference","text":"Modules = [PowerFlows]\nPublic = true\nPrivate = false","category":"page"},{"location":"api/public/#PowerFlows.PowerFlowData","page":"Public API Reference","title":"PowerFlows.PowerFlowData","text":"Structure containing all the data required for the evaluation of the power flows and angles, as well as these ones.\n\nArguments:\n\nbus_lookup::Dict{Int, Int}:       dictionary linking the system's bus number with the rows of either       \"powernetworkmatrix\" or \"auxnetworkmatrix\".\nbranch_lookup::Dict{String, Int}:       dictionary linking the branch name with the column name of either the       \"powernetworkmatrix\" or \"auxnetworkmatrix\".\nbus_activepower_injection::Matrix{Float64}:       \"(b, t)\" matrix containing the bus active power injection. b: number of       buses, t: number of time period.\nbus_reactivepower_injection::Matrix{Float64}:       \"(b, t)\" matrix containing the bus reactive power injection. b: number       of buses, t: number of time period.\nbus_activepower_withdrawals::Matrix{Float64}:       \"(b, t)\" matrix containing the bus reactive power withdrawals. b:       number of buses, t: number of time period.\nbus_reactivepower_withdrawals::Matrix{Float64}:       \"(b, t)\" matrix containing the bus reactive power withdrawals. b:       number of buses, t: number of time period.\nbus_reactivepower_bounds::Vector{Float64}: Upper and Lower bounds for the reactive supply       at each bus.\nbus_type::Vector{PSY.ACBusTypes}:       vector containing type of buses present in the system, ordered       according to \"bus_lookup\".\nbus_magnitude::Matrix{Float64}:       \"(b, t)\" matrix containing the bus magnitudes, ordered according to       \"bus_lookup\". b: number of buses, t: number of time period.\nbus_angles::Matrix{Float64}:       \"(b, t)\" matrix containing the bus angles, ordered according to       \"bus_lookup\". b: number of buses, t: number of time period.\nbranch_flow_values::Matrix{Float64}:       \"(br, t)\" matrix containing the power flows, ordered according to       \"branch_lookup\". br: number of branches, t: number of time period.\ntimestep_map::Dict{Int, S}:       dictonary mapping the number of the time periods (corresponding to the       column number of the previosly mentioned matrices) and their names.\nvalid_ix::Vector{Int}:       vector containing the indeces of not slack buses\npower_network_matrix::M:       matrix used for the evaluation of either the power flows or bus angles,       depending on the method considered.\naux_network_matrix::N:       matrix used for the evaluation of either the power flows or bus angles,       depending on the method considered.\nneighbors::Vector{Set{Int}}: Vector with the sets of adjacent buses.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerFlows.PowerFlowData-Tuple{ACPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.PowerFlowData","text":"Function for the definition of the PowerFlowData strucure given the System data, number of time periods to consider and their names. Calling this function will not evaluate the power flows and angles. NOTE: use it for AC power flow computations.\n\nArguments:\n\n::ACPowerFlow:       use ACPowerFlow() to evaluate the AC OPF.\nsys::PSY.System:       container storing the system data to consider in the PowerFlowData       structure.\ntime_steps::Int:       number of time periods to consider in the PowerFlowData structure. It       defines the number of columns of the matrices used to store data.       Default value = 1.\ntimestep_names::Vector{String}:       names of the time periods defines by the argmunet \"time_steps\". Default       value = String[].\ncheck_connectivity::Bool:       Perform connectivity check on the network matrix. Default value = true.\n\nWARNING: functions for the evaluation of the multi-period AC OPF still to be implemented.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.PowerFlowData-Tuple{DCPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.PowerFlowData","text":"Function for the definition of the PowerFlowData strucure given the System data, number of time periods to consider and their names. Calling this function will not evaluate the power flows and angles. NOTE: use it for DC power flow computations.\n\nArguments:\n\n::DCPowerFlow:       use DCPowerFlow() to store the ABA matrix as powernetworkmatrix and       the BA matrix as auxnetworkmatrix.\nsys::PSY.System:       container storing the system data to consider in the PowerFlowData       structure.\ntime_steps::Int:       number of time periods to consider in the PowerFlowData structure. It       defines the number of columns of the matrices used to store data.       Default value = 1.\ntimestep_names::Vector{String}:       names of the time periods defines by the argmunet \"time_steps\". Default       value = String[].\ncheck_connectivity::Bool:       Perform connectivity check on the network matrix. Default value = true.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.PowerFlowData-Tuple{PTDFDCPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.PowerFlowData","text":"Function for the definition of the PowerFlowData strucure given the System data, number of time periods to consider and their names. Calling this function will not evaluate the power flows and angles. NOTE: use it for DC power flow computations.\n\nArguments:\n\n::PTDFDCPowerFlow:       use PTDFDCPowerFlow() to store the PTDF matrix as powernetworkmatrix       and the ABA matrix as auxnetworkmatrix.\nsys::PSY.System:       container storing the system data to consider in the PowerFlowData       structure.\ntime_steps::Int:       number of time periods to consider in the PowerFlowData structure. It       defines the number of columns of the matrices used to store data.       Default value = 1.\ntimestep_names::Vector{String}:       names of the time periods defines by the argmunet \"time_steps\". Default       value = String[].\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.PowerFlowData-Tuple{vPTDFDCPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.PowerFlowData","text":"Function for the definition of the PowerFlowData strucure given the System data, number of time periods to consider and their names. Calling this function will not evaluate the power flows and angles. NOTE: use it for DC power flow computations.\n\nArguments:\n\n::PTDFDCPowerFlow:       use vPTDFDCPowerFlow() to store the Virtual PTDF matrix as       powernetworkmatrix and the ABA matrix as auxnetworkmatrix.\nsys::PSY.System:       container storing the system data to consider in the PowerFlowData       structure.\ntime_steps::Int:       number of time periods to consider in the PowerFlowData structure. It       defines the number of columns of the matrices used to store data.       Default value = 1.\ntimestep_names::Vector{String}:       names of the time periods defines by the argmunet \"time_steps\". Default       value = String[].\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_ac_powerflow!-Tuple{System}","page":"Public API Reference","title":"PowerFlows.solve_ac_powerflow!","text":"Solves a the power flow into the system and writes the solution into the relevant structs. Updates generators active and reactive power setpoints and branches active and reactive power flows (calculated in the From - To direction) (see flow_val)\n\nSupports passing NLsolve kwargs in the args. By default shows the solver trace.\n\nArguments available for nlsolve:\n\nget_connectivity::Bool: Checks if the network is connected. Default true\nmethod : See NLSolve.jl documentation for available solvers\nxtol: norm difference in x between two successive iterates under which convergence is declared. Default: 0.0.\nftol: infinite norm of residuals under which convergence is declared. Default: 1e-8.\niterations: maximum number of iterations. Default: 1_000.\nstore_trace: should a trace of the optimization algorithm's state be stored? Default: false.\nshow_trace: should a trace of the optimization algorithm's state be shown on STDOUT? Default: false.\nextended_trace: should additifonal algorithm internals be added to the state trace? Default: false.\n\nExamples\n\nsolve_ac_powerflow!(sys)\n# Passing NLsolve arguments\nsolve_ac_powerflow!(sys, method=:newton)\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_powerflow-Tuple{ACPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.solve_powerflow","text":"Similar to solve_powerflow!(sys) but does not update the system struct with results. Returns the results in a dictionary of dataframes.\n\nExamples\n\nres = solve_powerflow(sys)\n# Passing NLsolve arguments\nres = solve_powerflow(sys, method=:newton)\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_powerflow-Tuple{DCPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.solve_powerflow","text":"Evaluates the power flows on each system's branch by means of the ABA and BA matrices. Updates the PowerFlowData structure and returns a dictionary containing a DataFrame for the single timestep considered. The DataFrame containts the flows and angles related to the information stored in the PSY.System considered as input.\n\nArguments:\n\n::DCPowerFlow:       use DCPowerFlow() to evaluate the power flows according to the method       based on the ABA and BA matrices\nsys::PSY.System:       container gathering the system data used for the evaluation of flows       and angles.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_powerflow-Tuple{PTDFDCPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.solve_powerflow","text":"Evaluates the power flows on each system's branch by means of the PTDF matrix. Updates the PowerFlowData structure and returns a dictionary containing a DataFrame for the single timestep considered. The DataFrame containts the flows and angles related to the information stored in the PSY.System considered as input.\n\nArguments:\n\n::PTDFDCPowerFlow:       use PTDFDCPowerFlow() to evaluate the power flows according to the       method based on the PTDF matrix\nsys::PSY.System:       container gathering the system data used for the evaluation of flows       and angles.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_powerflow-Tuple{PowerFlowData, System}","page":"Public API Reference","title":"PowerFlows.solve_powerflow","text":"Evaluates the power flows on each system's branch by means of Virtual PTDF matrices. Updates the PowerFlowData structure \"data\" and returns a dictionary containing a number of DataFrames equal to the numeber of timestep considered in \"data\". Each DataFrame containts the flows and angles.\n\nArguments:\n\ndata::PTDFPowerFlowData:       PowerFlowData structure containing the system data per each timestep       considered, as well as the Virtual PTDF matrix.\nsys::PSY.System:       container gathering the system data.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_powerflow-Tuple{PowerFlowData{PowerNetworkMatrices.ABA_Matrix{Tuple{Vector{Int64}, Vector{Int64}}, Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}}, KLU.KLUFactorization{Float64, Int64}}, PowerNetworkMatrices.BA_Matrix{Tuple{Vector{Int64}, Vector{String}}, Tuple{Dict{Int64, Int64}, Dict{String, Int64}}}}, System}","page":"Public API Reference","title":"PowerFlows.solve_powerflow","text":"Evaluates the power flows on each system's branch by means of the ABA and BA matrices. Updates the PowerFlowData structure \"data\" and returns a dictionary containing a number of DataFrames equal to the numeber of timestep considered in \"data\". Each DataFrame containts the flows and angles.\n\nArguments:\n\ndata::ABAPowerFlowData:       PowerFlowData structure containing the system's data per each timestep       considered, as well as the ABA and BA matrices.\nsys::PSY.System:       container gathering the system data.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_powerflow-Tuple{PowerFlowData{PowerNetworkMatrices.PTDF{Tuple{Vector{Int64}, Vector{String}}, Tuple{Dict{Int64, Int64}, Dict{String, Int64}}, Matrix{Float64}}, PowerNetworkMatrices.ABA_Matrix{Tuple{Vector{Int64}, Vector{Int64}}, Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}}, KLU.KLUFactorization{Float64, Int64}}}, System}","page":"Public API Reference","title":"PowerFlows.solve_powerflow","text":"Evaluates the power flows on each system's branch by means of the PTDF matrix. Updates the PowerFlowData structure \"data\" and returns a dictionary containing a number of DataFrames equal to the numeber of timestep considered in \"data\". Each DataFrame containts the flows and angles.\n\nArguments:\n\ndata::PTDFPowerFlowData:       PowerFlowData structure containing the system's data per each timestep       considered, as well as the PTDF matrix.\nsys::PSY.System:       container gathering the system data.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.solve_powerflow-Tuple{vPTDFDCPowerFlow, System}","page":"Public API Reference","title":"PowerFlows.solve_powerflow","text":"Evaluates the power flows on each system's branch by means of the Virtual PTDF matrix. Updates the PowerFlowData structure \"data\" and returns a dictionary containing a number of DataFrames equal to the numeber of timestep considered in \"data\". The DataFrame containts the flows and angles related to the information stored in the PSY.System considered as input.\n\nArguments:\n\n::vPTDFDCPowerFlow:       use vPTDFDCPowerFlow() to evaluate the power flows according to the       method based on the Virtual PTDF matrix\nsys::PSY.System:       container gathering the system data used for the evaluation of flows       and angles.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.write_results-Tuple{ACPowerFlow, System, PowerFlowData{PowerNetworkMatrices.Ybus{Tuple{Vector{Int64}, Vector{Int64}}, Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}}}, Nothing}, Vector{Float64}}","page":"Public API Reference","title":"PowerFlows.write_results","text":"Returns a dictionary containing the AC power flow results.\n\nOnly single-period evaluation is supported at the moment for AC Power flows. Resulting dictionary will therefore feature just one key linked to one DataFrame.\n\nArguments:\n\n::ACPowerFlow:       use ACPowerFlow() storing AC power flow results.\nsys::PSY.System:       container storing the systam information.\nresult::Vector{Float64}:       vector containing the reults for one single time-period.\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerFlows.write_results-Tuple{PowerFlowData, System}","page":"Public API Reference","title":"PowerFlows.write_results","text":"Returns a dictionary containing the DC power flow results. Each key conresponds to the name of the considered time periods, storing a DataFrame with the OPF results.\n\nArguments:\n\ndata::Union{PTDFPowerFlowData, vPTDFPowerFlowData, ABAPowerFlowData}:       PowerFlowData strcuture containing power flows and bus angles.\nsys::PSY.System:       container storing the systam information.\n\n\n\n\n\n","category":"method"}]
}
