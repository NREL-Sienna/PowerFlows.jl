"""A wrapper for the factorization object for the Hessian matrix in the robust homotopy method.
Subtypes must implement: `symbolic_factor!`, `solve!`, and `modify_and_numeric_factor!`. 
The first two follow the same interface as `LinearSolverCache`. The last one modifies the 
Hessian to be invertible, then does a numeric factorization."""

abstract type HessianSolver end
