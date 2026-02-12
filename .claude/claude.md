# PowerFlows.jl Repository Guide

> **Development Guidelines:** For code development preferences, style conventions, and best practices for projects using Sienna, refer to [Sienna.md](./Sienna.md).

## Overview

PowerFlows.jl is a Julia package that provides an interface to multiple solution methods for the Power Flow problem, along with utilities commonly found in commercial software like Siemens PSS/e and GE PSLF. This package is part of the NREL-Sienna ecosystem focused on scalable power systems modeling and simulation.

## Main Objectives

1. **Multiple Solution Methods**: Provide a unified interface to various power flow solution algorithms:
   - DC Power Flow (linear approximation)
   - AC Power Flow using Newton-Raphson
   - AC Power Flow using Trust Region methods
   - AC Power Flow using Levenberg-Marquardt
   - Robust Homotopy methods for difficult-to-converge systems
   - PTDF-based DC Power Flow methods

2. **Commercial Software Compatibility**: Export power flow results in formats compatible with commercial software PSS/e. PSLF not supported yet.

3. **Integration with Sienna**: Seamless integration with PowerSystems.jl and PowerSimulations.jl (PSI) for comprehensive power systems analysis

4. **Performance**: Optimized linear algebra operations using sparse matrices and specialized solvers (KLU)

5. **Flexibility**: Support for various power system features including:
   - Multi-period power flow analysis
   - HVDC systems and LCC converters
   - Distributed slack bus models
   - Arc types and network reductions
   - Loss factor calculations
   - Stabiltiy Factor Calculations

## Repository Structure

### Core Source Code (`src/`)

#### Main Module
- **`PowerFlows.jl`**: Main module file containing exports and include statements

#### Power Flow Data and Setup
- **`PowerFlowData.jl`**: Core data structure holding power flow problem state
- **`power_flow_types.jl`**: Type definitions for different power flow methods
- **`initialize_power_flow_data.jl`**: Initialize power flow data structures from PowerSystems
- **`power_flow_setup.jl`**: Setup utilities for power flow problems
- **`power_flow_method.jl`**: Entry point for solving power flow

#### DC Power Flow
- **`solve_dc_power_flow.jl`**: DC power flow implementation and PTDF methods

#### AC Power Flow
- **`solve_ac_power_flow.jl`**: AC power flow implementations (Newton-Raphson, Trust Region)
- **`ac_power_flow_residual.jl`**: Residual (mismatch) calculations for AC power flow
- **`ac_power_flow_jacobian.jl`**: Jacobian matrix construction for Newton methods
- **`levenberg-marquardt.jl`**: Levenberg-Marquardt algorithm implementation

#### Robust Homotopy (`RobustHomotopy/`)
- **`robust_homotopy_method.jl`**: Robust homotopy continuation method for difficult cases
- **`homotopy_hessian.jl`**: Hessian calculations for homotopy methods
- **`HessianSolver/`**: Specialized solvers for Hessian systems (KLU, CHOLMOD, Cholesky)

#### Linear Algebra (`LinearSolverCache/`)
- **`linear_solver_cache.jl`**: Caching infrastructure for linear solvers
- **`klu_linear_solver.jl`**: KLU sparse direct solver interface

#### Utilities
- **`common.jl`**: Shared utility functions
- **`definitions.jl`**: Constants and type definitions
- **`state_indexing_helpers.jl`**: Helpers for indexing state variables
- **`lcc_parameters.jl`** & **`lcc_utils.jl`**: Line-commutated converter (HVDC) utilities
- **`powersystems_utils.jl`**: Integration utilities with PowerSystems.jl
- **`psi_utils.jl`**: Integration utilities with PowerSimulations.jl
- **`post_processing.jl`**: Result processing and calculations
- **`psse_export.jl`**: PSS/e format export functionality

### Tests (`test/`)

Comprehensive test suite organized by functionality:
- **`runtests.jl`**: Main test entry point
- **DC Power Flow**: `test_dc_power_flow.jl`, `test_multiperiod_dc_power_flow.jl`, `test_reduced_dc_power_flow.jl`
- **AC Power Flow**: `test_solve_power_flow.jl`, `test_multiperiod_ac_power_flow.jl`, `test_reduced_ac_power_flow.jl`
- **Specialized Features**:
  - `test_hvdc.jl`: HVDC systems
  - `test_distributed_slack.jl`: Distributed slack bus
  - `test_robust_power_flow.jl`: Robust homotopy methods
  - `test_iterative_methods.jl`: Iterative solution methods
- **Components**: `test_jacobian.jl`, `test_homotopy_hessian.jl`, `test_klu_linear_solver_cache.jl`
- **Utilities**: `test_psse_export.jl`, `test_post_processing.jl`, `test_power_flow_data.jl`, `test_loss_factors.jl`
- **Integration**: `test_psi_utils.jl`
- **Performance**: `performance/performance_test.jl`
- **Test Data**: `test_data/` contains validation datasets and reference results
- **Test Utilities**: `test_utils/` contains helper functions for testing

### Documentation (`docs/`)

- **`src/`**: Documentation source files (Markdown)
- **`build/`**: Generated documentation (not tracked in git)
- **`make.jl`**: Documentation build script
- Structure follows Diataxis framework:
  - **`tutorials/`**: Learning-oriented guides
  - **`how-tos/`**: Task-oriented guides
  - **`reference/`**: Information-oriented technical reference
    - `api/`: API documentation (public and internal)
    - `developers/`: Developer documentation
  - **`explanation/`**: Understanding-oriented explanations

### Scripts (`scripts/`)

- **`formatter/`**: Code formatting utilities using JuliaFormatter

### Configuration Files

- **`Project.toml`**: Package dependencies and metadata
- **`Manifest.toml`**: Exact versions of all dependencies (locked)
- **`codecov.yml`**: Code coverage configuration
- **`CONTRIBUTING.md`**: Contribution guidelines
- **`LICENSE`**: BSD license

## Key Concepts and Workflows

### Power Flow Problem

The power flow (load flow) problem solves for steady-state voltages (magnitude and angle) at all buses in a power network given:
- Generation levels
- Load demands
- Network topology and parameters

### Typical Usage Pattern

1. **Create System**: Use PowerSystems.jl to create or load a `System` object
2. **Choose Method**: Select appropriate power flow method (DC, AC with specific algorithm)
3. **Solve**: Call `solve_power_flow` or `solve_power_flow!` (in-place)
4. **Analyze Results**: Extract voltages, flows, losses, etc. from the results
5. **Export** (optional): Export results in PSS/e format using `PSSEExporter`

### Solution Methods Selection Guide

- **DCPowerFlow**: Fast linear approximation, good for screening studies
- **NewtonRaphsonACPowerFlow**: Standard method, fast convergence for well-conditioned systems
- **TrustRegionACPowerFlow**: More robust than Newton-Raphson, handles ill-conditioned cases better
- **LevenbergMarquardtACPowerFlow**: Robust nonlinear solver, good for difficult cases
- **RobustHomotopyPowerFlow**: Most robust method for hard-to-converge or non-convergent cases
- **PTDFDCPowerFlow** / **vPTDFDCPowerFlow**: DC power flow with pre-computed power transfer distribution factors

### Integration with PowerSimulations.jl (PSI)

PowerFlows is used within PowerSimulations.jl for:
- Initialization of dyanmic simulation problems
- Network model validation
- Post-processing of optimal power flow results

The `psi_utils.jl` file provides the integration layer.

## Development Workflow

1. **Setup**: Clone repository, instantiate environment with `]instantiate` in Julia REPL
2. **Code**: Follow Sienna conventions (see [Sienna.md](./Sienna.md))
3. **Format**: Run formatter from `scripts/formatter/` before committing
4. **Test**: Run tests with `]test PowerFlows` or specific test files
5. **Document**: Add docstrings following Diataxis framework
6. **Submit**: Create pull request following CONTRIBUTING.md

## Active Development Areas

Based on the current PR context (Refactor: Updates on PSS/e Exporter, testing DCPF validation on EI/WECC):
- **PSS/e Export Functionality**: Ongoing refactoring and enhancement of power system export capabilities
- **DC Power Flow Validation**: Testing and validation on large-scale systems (Eastern Interconnect, WECC)
- **Exporter Architecture**: Improving the design and maintainability of export functionality

## Dependencies

### Key External Packages
- **PowerSystems.jl**: Power system data models and network representation
- **PowerNetworkMatrices.jl**: Network matrix formulations (admittance, incidence, PTDF)
- **InfrastructureSystems.jl**: Core infrastructure shared across Sienna packages
- **KLU.jl**: Sparse direct solver (from SuiteSparse)
- **DataFrames.jl**: Tabular data handling for results

## Resources

- **Documentation**: https://nrel-sienna.github.io/PowerFlows.jl/dev/
- **Issues**: GitHub Issues for bug reports and feature requests
- **Style Guide**: https://nrel-sienna.github.io/InfrastructureSystems.jl/stable/style/

## Version Information

- **Julia Compatibility**: 1.10+
- **Status**: Active development (API subject to change)

---

*This document is maintained to help developers (human and AI) understand and contribute to PowerFlows.jl effectively.*
