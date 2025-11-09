# CI/CD Workflows

This directory contains the GitHub Actions workflows for PowerFlows.jl.

## Workflows

### `main-tests.yml`
Runs on every push to the `main` branch and on a scheduled basis.

**Jobs:**
- `test`: Standard test suite across multiple OS (Ubuntu, Windows, macOS) and Julia versions
- `test-cuda`: CUDA-enabled tests on Ubuntu with CUDA toolkit installed

### `pr_testing.yml`
Runs on pull requests (opened, synchronized, or reopened).

**Jobs:**
- `test`: Standard test suite across multiple OS (Ubuntu, Windows, macOS)
- `test-cuda`: CUDA-enabled tests on Ubuntu with CUDA toolkit installed

### `docs.yml`
Builds and deploys documentation.

### `format-check.yml`
Checks code formatting.

### `TagBot.yml`
Automated tagging and release creation.

## CUDA Testing

The `test-cuda` job in both main and PR testing workflows provides GPU-accelerated linear solver testing:

**Features:**
- Installs CUDA Toolkit 12.2.0 on Ubuntu runners
- Installs CUDA.jl for Julia GPU computing
- Runs all tests with CUDA extension loaded
- Provides separate code coverage reporting for CUDA code paths

**Note:** CUDA tests use CPU-based CUDA simulation since GitHub Actions runners don't have physical GPUs. The tests verify:
- CUDA extension loads correctly
- CUSOLVERLinSolveCache can be instantiated
- Linear solver interface is properly implemented
- Integration with AC and DC power flow methods

**Environment Variables:**
- `JULIA_CUDA_USE_BINARYBUILDER="false"`: Uses system CUDA installation instead of BinaryBuilder artifacts

## Code Coverage

Code coverage is reported to Codecov with separate flags:
- `unittests`: Standard test coverage
- `unittests-cuda`: CUDA-specific test coverage

This allows tracking coverage for both CPU-only and GPU-accelerated code paths.
