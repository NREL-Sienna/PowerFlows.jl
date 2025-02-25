"""
    struct ACPowerFlowJacobian

A struct that represents the Jacobian matrix for AC power flow calculations.

This struct uses the functor pattern, meaning instances of `ACPowerFlowJacobian` store the data (Jacobian matrix) internally
and can be called as a function at the same time. Calling the instance as a function updates the stored Jacobian matrix.


# Fields
- `data::ACPowerFlowData`: The grid model data used for power flow calculations.
- `Jf!::Function`: A function that calculates the Jacobian matrix inplace.
- `Jv::SparseArrays.SparseMatrixCSC{Float64, Int32}`: The Jacobian matrix, which is updated by the function `Jf!`.
"""
struct ACPowerFlowJacobian
    data::ACPowerFlowData
    Jf!::Function   # This is the function that calculates the Jacobian matrix and updates Jv inplace
    Jv::SparseArrays.SparseMatrixCSC{Float64, Int32}  # This is the Jacobian matrix, that is updated by the function Jf
end

"""
    (J::ACPowerFlowJacobian)(time_step::Int64)

Update the Jacobian matrix `Jv` using the function `Jf!` and the provided data and time step.

Defining this method allows an instance of `ACPowerFlowJacobian` to be called as a function, following the functor pattern.

# Arguments
- `time_step::Int64`: The time step for the calculations.

# Example
```julia
J = ACPowerFlowJacobian(data, time_step)
J(time_step)  # Updates the Jacobian matrix Jv
```
"""
function (J::ACPowerFlowJacobian)(time_step::Int64)
    J.Jf!(J.Jv, J.data, time_step)
    return
end

"""
    (J::ACPowerFlowJacobian)(J::SparseArrays.SparseMatrixCSC{Float64, Int32}, time_step::Int64)

Use the `ACPowerFlowJacobian` to update the provided Jacobian matrix `J` inplace.

Update the internally stored Jacobian matrix `Jv` using the function `Jf!` and the provided data and time step, and write the updated Jacobian values to `J`.

This method allows an instance of ACPowerFlowJacobian to be called as a function, following the functor pattern.

# Arguments
- `J::SparseArrays.SparseMatrixCSC{Float64, Int32}``: A sparse matrix to be updated with new values of the Jacobian matrix.
- `time_step::Int64`: The time step for the calculations.

# Example
```julia
J = ACPowerFlowJacobian(data, time_step)
Jv = SparseArrays.sparse(Float64[], Int32[], Int32[])
J(Jv, time_step)  # Updates the Jacobian matrix Jv and writes it to J
```
"""
function (J::ACPowerFlowJacobian)(
    Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    time_step::Int64,
)
    J.Jf!(J.Jv, J.data, time_step)
    copyto!(Jv, J.Jv)
    return
end

"""
    ACPowerFlowJacobian(data::ACPowerFlowData, time_step::Int64) -> ACPowerFlowJacobian

This is the constructor for ACPowerFlowJacobian.
Create an `ACPowerFlowJacobian` instance. As soon as the instance is created, it already has the Jacobian matrix structure initialized and its values updated, stored internally as Jv.
The data instance is stored internally and used to update the Jacobian matrix because the structure of the Jacobian matrix is tied to the data.
Changing the data requires creating a new instance of `ACPowerFlowJacobian`.

# Arguments
- `data::ACPowerFlowData`: The data used for power flow calculations.
- `time_step::Int64`: The time step for the calculations.

# Returns
- `ACPowerFlowJacobian`: An instance of `ACPowerFlowJacobian`.

#Example
```julia
J = ACPowerFlowJacobian(data, time_step)  # Creates an instance J of ACPowerFlowJacobian, with the Jacobian matrix stored internally as J.Jv initialized and updated.
J(time_step)  # Updates the Jacobian matrix stored internally in J (J.Jv) with the latest state of the `data` (`ACPowerFlowData` instance) and the provided time step.
J.Jv  # Access the Jacobian matrix stored internally in J.
```
"""
function ACPowerFlowJacobian(data::ACPowerFlowData, time_step::Int64)
    # Create the initial Jacobian matrix structure - a sparse matrix with structural zeros that will be updated by the function Jf!
    # It has the same structure as the expected Jacobian matrix.
    Jv0 = _create_jacobian_matrix_structure(data, time_step)
    # _update_jacobian_matrix_values!(Jv0, data, time_step)  # We just initialize the structure here, evaluation must happen later
    return ACPowerFlowJacobian(data, _update_jacobian_matrix_values!, Jv0)
end

"""
Create the Jacobian matrix structure for a reference bus (REF). Ignoring this because we fill all four values even for PV buses with 
    structiural zeros using the same function as for PQ buses.
"""
function _create_jacobian_matrix_structure_bus!(rows::Vector{Int32},
    columns::Vector{Int32},
    values::Vector{Float64},
    bus_from::Int,
    bus_to::Int,
    row_from_p::Int,
    row_from_q::Int,
    col_to_vm::Int,
    col_to_va::Int,
    ::Val{PSY.ACBusTypes.REF})
    if bus_from == bus_to
        # Active PF w/r Local Active Power
        push!(rows, row_from_p)
        push!(columns, col_to_vm)
        push!(values, 0.0)
        # Rective PF w/r Local Reactive Power
        push!(rows, row_from_q)
        push!(columns, col_to_va)
        push!(values, 0.0)
    end
    return nothing
end

"""
Create the Jacobian matrix structure for a PV bus. Ignoring this because we fill all four values even for PV buses with 
    structiural zeros using the same function as for PQ buses.
"""
function _create_jacobian_matrix_structure_bus!(rows::Vector{Int32},
    columns::Vector{Int32},
    values::Vector{Float64},
    bus_from::Int,
    bus_to::Int,
    row_from_p::Int,
    row_from_q::Int,
    col_to_vm::Int,
    col_to_va::Int,
    ::Val{PSY.ACBusTypes.PV})
    # Active PF w/r Voltage Angle
    push!(rows, row_from_p)
    push!(columns, col_to_va)
    push!(values, 0.0)
    # Reactive PF w/r Voltage Angle
    push!(rows, row_from_q)
    push!(columns, col_to_va)
    push!(values, 0.0)
    if bus_from == bus_to
        # Reactive PF w/r Local Reactive Power
        push!(rows, row_from_q)
        push!(columns, col_to_vm)
        push!(values, 0.0)
    end
    return nothing
end

"""
Create the Jacobian matrix structure for a PQ bus. Using this for all buses because 
    a) for REF buses it doesn't matter if there are 2 values or 4 values - there are not many of them in the grid
    b) for PV buses we fill all four values because we can have a PV -> PQ transition and then we need to fill all four values
"""
function _create_jacobian_matrix_structure_bus!(rows::Vector{Int32},
    columns::Vector{Int32},
    values::Vector{Float64},
    bus_from::Int,
    bus_to::Int,
    row_from_p::Int,
    row_from_q::Int,
    col_to_vm::Int,
    col_to_va::Int,
    # ::Val{PSY.ACBusTypes.PQ}
    )
    # Active PF w/r Voltage Magnitude
    push!(rows, row_from_p)
    push!(columns, col_to_vm)
    push!(values, 0.0)
    # Reactive PF w/r Voltage Magnitude
    push!(rows, row_from_q)
    push!(columns, col_to_vm)
    push!(values, 0.0)
    # Active PF w/r Voltage Angle
    push!(rows, row_from_p)
    push!(columns, col_to_va)
    push!(values, 0.0)
    # Reactive PF w/r Voltage Angle
    push!(rows, row_from_q)
    push!(columns, col_to_va)
    push!(values, 0.0)
    return nothing
end

"""
    _create_jacobian_matrix_structure(data::ACPowerFlowData, time_step::Int64) -> SparseMatrixCSC{Float64, Int32}

Create the structure of the Jacobian matrix for an AC power flow problem. Inputs are the grid model as an instance of `ACPowerFlowData` at a given time step.

# Arguments
- `data::ACPowerFlowData`: The power flow model.
- `time_step::Int64`: The specific time step for which the Jacobian matrix structure is created.

# Returns
- `SparseMatrixCSC{Float64, Int32}`: A sparse matrix with structural zeros representing the structure of the Jacobian matrix.

# Description
This function initializes the structure of the Jacobian matrix for an AC power flow problem. 
The Jacobian matrix is used in power flow analysis to represent the partial derivatives of bus active and reactive power injections with respect to bus voltage magnitudes and angles.

Unlike some commonly used approaches where the Jacobian matrix is constructed as four submatrices, each grouping values for the four types of partial derivatives,
this function groups the partial derivatives by bus.
The structure is organized as groups of 4 values per bus, corresponding to the following partial derivatives:

| ∂P₁/∂V₁  | ∂P₁/∂θ₁ | ...      | ...     | ...     | ...     | ...     |
| ∂Q₁/∂V₁  | ∂Q₁/∂θ₁ | ...      | ...     | ...     | ...     | ...     | 
| ...      |         | ∂P₂/∂V₂  | ∂P₂/∂θ₂ | ...     | ...     | ...     |
| ...      |         | ∂Q₂/∂V₂  | ∂Q₂/∂θ₂ | ...     | ...     | ...     |
| ...      | ...     | ...      | ...     | ...     | ∂Pₙ/∂Vₙ  | ∂Pₙ/∂θₙ |
| ...      | ...     | ...      | ...     | ...     | ∂Qₙ/∂Vₙ  | ∂Qₙ/∂θₙ |


This approach is more memory-efficient. Furthermore, this structure results in a more efficient factorization because the values are more likely to be grouped close to the diagonal.
Refer to Electric Energy Systems: Analysis and Operation by Antonio Gomez-Exposito and Fernando L. Alvarado for more details.

The function initializes three arrays (`rows`, `columns`, and `values`) to store the row indices, column indices, and values of the non-zero elements of the Jacobian matrix, respectively.

For each bus in the system, the function iterates over its neighboring buses and determines the type of each neighboring bus (`REF`, `PV`, or `PQ`). 
Depending on the bus type, the function adds the appropriate entries to the Jacobian matrix structure.

- For `REF` buses, entries are added for local active and reactive power.
- For `PV` buses, entries are added for active and reactive power with respect to angle, and for local reactive power.
- For `PQ` buses, entries are added for active and reactive power with respect to voltage magnitude and angle.

Finally, the function constructs a sparse matrix from the collected indices and values and returns it.
"""
function _create_jacobian_matrix_structure(data::ACPowerFlowData, time_step::Int64)
    # Create Jacobian structure
    # Initialize arrays to store the row indices, column indices, and values of the non-zero elements of the Jacobian matrix
    rows = Int32[]      # I
    columns = Int32[]   # J
    values = Float64[]  # V

    num_buses = first(size(data.bus_type))

    for bus_from in 1:num_buses
        row_from_p = 2 * bus_from - 1  # Row index for the value that is related to active power
        row_from_q = 2 * bus_from      # Row index for the value that is related to reactive power
        for bus_to in data.neighbors[bus_from]
            col_to_vm = 2 * bus_to - 1  # Column index for the value related to voltage magnitude
            col_to_va = 2 * bus_to      # Column index for the value related to voltage angle
            # We ignore the bus type and initialize the structure as if all buses were PQ - 
            # mainly because we can have a PV -> PQ transition, and the number of REF buses is small
            # bus_type = data.bus_type[bus_to, time_step]
            _create_jacobian_matrix_structure_bus!(
                rows,
                columns,
                values,
                bus_from,
                bus_to,
                row_from_p,
                row_from_q,
                col_to_vm,
                col_to_va,
                # Val(bus_type),
            )
        end
    end
    Jv0 = SparseArrays.sparse(rows, columns, values)
    return Jv0
end

function _set_entries_for_neighbor(Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Yb::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    Vm::Vector{Float64},
    θ::Vector{Float64},
    bus_from::Int, bus_to::Int,
    row_from_p::Int, row_from_q::Int,
    col_to_vm::Int, col_to_va::Int,
    bus_from_neighbors::Set{Int},
    ::Val{PSY.ACBusTypes.REF})
    # State variables are Active and Reactive Power Generated
    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - x[2*i-1]
    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
    # x does not appear in p_flow and q_flow
    if bus_from == bus_to
        Jv[row_from_p, col_to_vm] = -1.0
        Jv[row_from_q, col_to_va] = -1.0
    end
end

function _set_entries_for_neighbor(Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Yb::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    Vm::Vector{Float64},
    θ::Vector{Float64},
    bus_from::Int, bus_to::Int,
    row_from_p::Int, row_from_q::Int,
    col_to_vm::Int, col_to_va::Int,
    bus_from_neighbors::Set{Int},
    ::Val{PSY.ACBusTypes.PV})
    # State variables are Reactive Power Generated and Voltage Angle
    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - p_gen[i]
    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
    # x[2*i] (associated with q_gen) does not appear in q_flow
    # x[2*i] (associated with q_gen) does not appear in the active power balance
    if bus_from == bus_to
        # Jac: Reactive PF against local active power
        Jv[row_from_q, col_to_vm] = -1.0
        # Jac: Active PF against same Angle: θ[bus_from] =  θ[bus_to]
        Jv[row_from_p, col_to_va] =
            Vm[bus_from] * sum(
                Vm[k] * (
                    real(Yb[bus_from, k]) * -sin(θ[bus_from] - θ[k]) +
                    imag(Yb[bus_from, k]) * cos(θ[bus_from] - θ[k])
                ) for k in bus_from_neighbors if k != bus_from
            )
        # Jac: Reactive PF against same Angle: θ[bus_from] = θ[bus_to]
        Jv[row_from_q, col_to_va] =
            Vm[bus_from] * sum(
                Vm[k] * (
                    real(Yb[bus_from, k]) * cos(θ[bus_from] - θ[k]) -
                    imag(Yb[bus_from, k]) * -sin(θ[bus_from] - θ[k])
                ) for k in bus_from_neighbors if k != bus_from
            )
    else
        g_ij = real(Yb[bus_from, bus_to])
        b_ij = imag(Yb[bus_from, bus_to])
        # Jac: Active PF against other angles θ[bus_to]
        Jv[row_from_p, col_to_va] =
            Vm[bus_from] *
            Vm[bus_to] *
            (g_ij * sin(θ[bus_from] - θ[bus_to]) + b_ij * -cos(θ[bus_from] - θ[bus_to]))
        # Jac: Reactive PF against other angles θ[bus_to]
        Jv[row_from_q, col_to_va] =
            Vm[bus_from] *
            Vm[bus_to] *
            (g_ij * -cos(θ[bus_from] - θ[bus_to]) - b_ij * sin(θ[bus_from] - θ[bus_to]))
    end
end

function _set_entries_for_neighbor(Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Yb::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    Vm::Vector{Float64},
    θ::Vector{Float64},
    bus_from::Int, bus_to::Int,
    row_from_p::Int, row_from_q::Int,
    col_to_vm::Int, col_to_va::Int,
    bus_neighbors::Set{Int},
    ::Val{PSY.ACBusTypes.PQ})
    # State variables are Voltage Magnitude and Voltage Angle
    # Everything appears in everything
    if bus_from == bus_to
        # Jac: Active PF against same voltage magnitude Vm[bus_from]
        Jv[row_from_p, col_to_vm] =
            2 * real(Yb[bus_from, bus_to]) * Vm[bus_from] + sum(
                Vm[k] * (
                    real(Yb[bus_from, k]) * cos(θ[bus_from] - θ[k]) +
                    imag(Yb[bus_from, k]) * sin(θ[bus_from] - θ[k])
                ) for k in bus_neighbors if k != bus_from
            )
        # Jac: Active PF against same angle θ[bus_from]
        Jv[row_from_p, col_to_va] =
            Vm[bus_from] * sum(
                Vm[k] * (
                    real(Yb[bus_from, k]) * -sin(θ[bus_from] - θ[k]) +
                    imag(Yb[bus_from, k]) * cos(θ[bus_from] - θ[k])
                ) for k in bus_neighbors if k != bus_from
            )

        # Jac: Reactive PF against same voltage magnitude Vm[bus_from]
        Jv[row_from_q, col_to_vm] =
            -2 * imag(Yb[bus_from, bus_to]) * Vm[bus_from] + sum(
                Vm[k] * (
                    real(Yb[bus_from, k]) * sin(θ[bus_from] - θ[k]) -
                    imag(Yb[bus_from, k]) * cos(θ[bus_from] - θ[k])
                ) for k in bus_neighbors if k != bus_from
            )
        # Jac: Reactive PF against same angle θ[bus_from]
        Jv[row_from_q, col_to_va] =
            Vm[bus_from] * sum(
                Vm[k] * (
                    real(Yb[bus_from, k]) * cos(θ[bus_from] - θ[k]) -
                    imag(Yb[bus_from, k]) * -sin(θ[bus_from] - θ[k])
                ) for k in bus_neighbors if k != bus_from
            )
    else
        g_ij = real(Yb[bus_from, bus_to])
        b_ij = imag(Yb[bus_from, bus_to])
        # Jac: Active PF w/r to different voltage magnitude Vm[bus_to]
        Jv[row_from_p, col_to_vm] =
            Vm[bus_from] *
            (g_ij * cos(θ[bus_from] - θ[bus_to]) + b_ij * sin(θ[bus_from] - θ[bus_to]))
        # Jac: Active PF w/r to different angle θ[bus_to]
        Jv[row_from_p, col_to_va] =
            Vm[bus_from] *
            Vm[bus_to] *
            (g_ij * sin(θ[bus_from] - θ[bus_to]) + b_ij * -cos(θ[bus_from] - θ[bus_to]))

        # Jac: Reactive PF w/r to different voltage magnitude Vm[bus_to]
        Jv[row_from_q, col_to_vm] =
            Vm[bus_from] *
            (g_ij * sin(θ[bus_from] - θ[bus_to]) - b_ij * cos(θ[bus_from] - θ[bus_to]))
        # Jac: Reactive PF w/r to different angle θ[bus_to]
        Jv[row_from_q, col_to_va] =
            Vm[bus_from] *
            Vm[bus_to] *
            (g_ij * -cos(θ[bus_from] - θ[bus_to]) - b_ij * sin(θ[bus_from] - θ[bus_to]))
    end
end

"""Used to update Jv based on the bus voltages, angles, etc. in data."""
function _update_jacobian_matrix_values!(
    Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    data::ACPowerFlowData,
    time_step::Int64,
)
    Yb = data.power_network_matrix.data
    Vm = data.bus_magnitude[:, time_step]
    θ = data.bus_angles[:, time_step]
    num_buses = first(size(data.bus_type))
    for bus_from in 1:num_buses
        row_from_p = 2 * bus_from - 1
        row_from_q = 2 * bus_from

        for bus_to in data.neighbors[bus_from]
            col_to_vm = 2 * bus_to - 1
            col_to_va = 2 * bus_to
            bus_type = data.bus_type[bus_to, time_step]
            _set_entries_for_neighbor(Jv, Yb, Vm, θ,
                bus_from, bus_to,
                row_from_p, row_from_q,
                col_to_vm, col_to_va,
                data.neighbors[bus_from], Val(bus_type))
        end
    end
    return
end
