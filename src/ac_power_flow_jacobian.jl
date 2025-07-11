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
    # Create the initial Jacobian matrix structure - a sparse matrix with structural zeros
    # that will be updated by the function Jf! It has the same structure as the expected
    # Jacobian matrix.
    Jv0 = _create_jacobian_matrix_structure(data, time_step)
    # We just initialize the structure here, evaluation must happen later
    return ACPowerFlowJacobian(data, _update_jacobian_matrix_values!, Jv0)
end

"""
Create the Jacobian matrix structure for a reference bus (REF). Ignoring this because we fill all four values even for PV buses with 
    structural zeros using the same function as for PQ buses.
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
        # Reactive PF w/r Local Reactive Power
        push!(rows, row_from_q)
        push!(columns, col_to_va)
        push!(values, 0.0)
    end
    return nothing
end

"""
Create the Jacobian matrix structure for a PV bus. Ignoring this because we fill all four values even for PV buses with 
    structural zeros using the same function as for PQ buses.
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
this function groups the partial derivatives by bus. The structure is organized as groups of 4 values per bus. See the example below for details.

This approach is more memory-efficient. Furthermore, this structure results in a more efficient factorization because the values are more likely to be grouped close to the diagonal.
Refer to Electric Energy Systems: Analysis and Operation by Antonio Gomez-Exposito and Fernando L. Alvarado for more details.

The function initializes three arrays (`rows`, `columns`, and `values`) to store the row indices, column indices, and values of the non-zero elements of the Jacobian matrix, respectively.

For each bus in the system, the function iterates over its neighboring buses and determines the type of each neighboring bus (`REF`, `PV`, or `PQ`). 
Depending on the bus type, the function adds the appropriate entries to the Jacobian matrix structure.

- For `REF` buses, entries are added for local active and reactive power.
- For `PV` buses, entries are added for active and reactive power with respect to angle, and for local reactive power.
- For `PQ` buses, entries are added for active and reactive power with respect to voltage magnitude and angle.

For example, suppose we have a system with 3 buses: bus 1 is `REF`, bus 2 is `PV`, and bus 3 is `PQ`.
Let ΔPⱼ, ΔQⱼ be the active, reactive power balance at the `j`th bus. Let Pⱼ and Qⱼ be the
active and reactive power generated at the `j`th bus (`REF` and `PV` only). Then the state vector is
[P₁, Q₁, Q₂, θ₂, V₃, θ₃], and the Jacobian matrix is

| ∂ΔP₁/∂P₁ | ∂ΔP₁/∂Q₁ | ∂ΔP₁/∂Q₂ | ∂ΔP₁/∂θ₂ | ∂ΔP₁/∂V₃ | ∂ΔP₁/∂θ₃ |  
| ∂ΔQ₁/∂P₁ | ∂ΔQ₁/∂Q₁ | ∂ΔQ₁/∂Q₂ | ∂ΔQ₁/∂θ₂ | ∂ΔQ₁/∂V₃ | ∂ΔQ₁/∂θ₃ |
| ∂ΔP₂/∂P₁ | ∂ΔP₂/∂Q₁ | ∂ΔP₂/∂Q₂ | ∂ΔP₂/∂θ₂ | ∂ΔP₂/∂V₃ | ∂ΔP₂/∂θ₃ |
| ∂ΔQ₂/∂P₁ | ∂ΔQ₂/∂Q₁ | ∂ΔQ₂/∂Q₂ | ∂ΔQ₂/∂θ₂ | ∂ΔQ₂/∂V₃ | ∂ΔQ₂/∂θ₃ |
| ∂ΔP₃/∂P₁ | ∂ΔP₃/∂Q₁ | ∂ΔP₃/∂Q₂ | ∂ΔP₃/∂θ₂ | ∂ΔP₃/∂V₃ | ∂ΔP₃/∂θ₃ |
| ∂ΔQ₃/∂P₁ | ∂ΔQ₃/∂Q₁ | ∂ΔQ₃/∂Q₂ | ∂ΔQ₃/∂θ₂ | ∂ΔQ₃/∂V₃ | ∂ΔQ₃/∂θ₃ |

In reality, for large networks, this matrix would be sparse, and each 4x4 block would only be nonzero
when there's a line between the respective buses.

Finally, the function constructs a sparse matrix from the collected indices and values and returns it.
"""
function _create_jacobian_matrix_structure(data::ACPowerFlowData, time_step::Int64)
    # Create Jacobian structure
    # Initialize arrays to store the row indices, column indices, and values of the non-zero elements of the Jacobian matrix
    rows = Int32[]      # I
    columns = Int32[]   # J
    values = Float64[]  # V

    num_buses = first(size(data.bus_type))

    num_lines = length(get_branch_lookup(data))
    sizehint!(rows, 4 * num_lines)
    sizehint!(columns, 4 * num_lines)
    sizehint!(values, 4 * num_lines)

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

function _set_entries_for_neighbor(::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Y_from_to::ComplexF32,
    Vm_from::Float64,
    Vm_to::Float64,
    θ_from_to::Float64,
    ::Int,
    ::Int,
    ::Int,
    ::Int,
    ∂P∂θ_from::Base.RefValue{Float64},
    ∂Q∂θ_from::Base.RefValue{Float64},
    ∂P∂V_from::Base.RefValue{Float64},
    ∂Q∂V_from::Base.RefValue{Float64},
    ::Val{PSY.ACBusTypes.REF})
    # State variables are Active and Reactive Power Generated
    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - x[2*i-1]
    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
    # x does not appear in p_flow and q_flow
    g_ij, b_ij = real(Y_from_to), imag(Y_from_to)
    # still need to do diagonal terms: those are based off
    # the bus type of from_bus, when we're dispatching on bustype of to_bus.
    ∂P∂θ_from[] -= Vm_from * Vm_to * (g_ij * sin(θ_from_to) - b_ij * cos(θ_from_to))
    ∂Q∂θ_from[] -= Vm_from * Vm_to * (-g_ij * cos(θ_from_to) - b_ij * sin(θ_from_to))
    ∂P∂V_from[] += Vm_to * (g_ij * cos(θ_from_to) + b_ij * sin(θ_from_to))
    ∂Q∂V_from[] += Vm_to * (g_ij * sin(θ_from_to) - b_ij * cos(θ_from_to))
    return
end

function _set_entries_for_neighbor(Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Y_from_to::ComplexF32,
    Vm_from::Float64,
    Vm_to::Float64,
    θ_from_to::Float64,
    row_from_p::Int,
    row_from_q::Int,
    ::Int,
    col_to_va::Int,
    ∂P∂θ_from::Base.RefValue{Float64},
    ∂Q∂θ_from::Base.RefValue{Float64},
    ∂P∂V_from::Base.RefValue{Float64},
    ∂Q∂V_from::Base.RefValue{Float64},
    ::Val{PSY.ACBusTypes.PV},
)
    # State variables are Reactive Power Generated and Voltage Angle
    # F[2*i-1] := p[i] = p_flow[i] + p_load[i] - p_gen[i]
    # F[2*i] := q[i] = q_flow[i] + q_load[i] - x[2*i]
    # x[2*i] (associated with q_gen) does not appear in q_flow
    # x[2*i] (associated with q_gen) does not appear in the active power balance
    g_ij, b_ij = real(Y_from_to), imag(Y_from_to)
    # Jac: Active PF against other angles θ[bus_to]
    p_va_common_term = Vm_from * Vm_to * (g_ij * sin(θ_from_to) - b_ij * cos(θ_from_to))
    Jv[row_from_p, col_to_va] = p_va_common_term
    ∂P∂θ_from[] -= p_va_common_term
    # Jac: Reactive PF w/r to different angle θ[bus_to]
    q_va_common_term = Vm_from * Vm_to * (-g_ij * cos(θ_from_to) - b_ij * sin(θ_from_to))
    Jv[row_from_q, col_to_va] = q_va_common_term
    ∂Q∂θ_from[] -= q_va_common_term

    # still need to do all diagonal terms: those are based off
    # the bus type of from_bus, when we're dispatching on bustype of to_bus.
    ∂P∂V_from[] += Vm_to * (g_ij * cos(θ_from_to) + b_ij * sin(θ_from_to))
    ∂Q∂V_from[] += Vm_to * (g_ij * sin(θ_from_to) - b_ij * cos(θ_from_to))
    return
end

function _set_entries_for_neighbor(Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    Y_from_to::ComplexF32,
    Vm_from::Float64,
    Vm_to::Float64,
    θ_from_to::Float64,
    row_from_p::Int,
    row_from_q::Int,
    col_to_vm::Int,
    col_to_va::Int,
    ∂P∂θ_from::Base.RefValue{Float64},
    ∂Q∂θ_from::Base.RefValue{Float64},
    ∂P∂V_from::Base.RefValue{Float64},
    ∂Q∂V_from::Base.RefValue{Float64},
    ::Val{PSY.ACBusTypes.PQ},
)
    # State variables are Voltage Magnitude and Voltage Angle
    # both state variables appear in both outputs.
    g_ij, b_ij = real(Y_from_to), imag(Y_from_to)
    # Active PF w/r to different voltage magnitude Vm[bus_to]
    p_vm_common_term = g_ij * cos(θ_from_to) + b_ij * sin(θ_from_to)
    Jv[row_from_p, col_to_vm] = Vm_from * p_vm_common_term
    ∂P∂V_from[] += Vm_to * p_vm_common_term
    # Active PF w/r to different angle θ[bus_to]
    p_va_common_term = Vm_from * Vm_to * (g_ij * sin(θ_from_to) - b_ij * cos(θ_from_to))
    Jv[row_from_p, col_to_va] = p_va_common_term
    ∂P∂θ_from[] -= p_va_common_term
    # Reactive PF w/r to different voltage magnitude Vm[bus_to]
    q_vm_common_term = g_ij * sin(θ_from_to) - b_ij * cos(θ_from_to)
    Jv[row_from_q, col_to_vm] = Vm_from * q_vm_common_term
    ∂Q∂V_from[] += Vm_to * q_vm_common_term
    # Jac: Reactive PF w/r to different angle θ[bus_to]
    q_va_common_term = Vm_from * Vm_to * (-g_ij * cos(θ_from_to) - b_ij * sin(θ_from_to))
    Jv[row_from_q, col_to_va] = q_va_common_term
    ∂Q∂θ_from[] -= q_va_common_term
    return
end

"""Used to update Jv based on the bus voltages, angles, etc. in data."""
function _update_jacobian_matrix_values!(
    Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    data::ACPowerFlowData,
    time_step::Int64,
)
    Yb = data.power_network_matrix.data
    Vm = view(data.bus_magnitude, :, time_step)
    θ = view(data.bus_angles, :, time_step)
    num_buses = first(size(data.bus_type))

    for bus_from in 1:num_buses
        row_from_p = 2 * bus_from - 1
        row_from_q = 2 * bus_from

        # the diagonal terms: e.g. ∂P_from/∂θ_from
        ∂P∂θ_from = Base.RefValue{Float64}(0.0)
        ∂Q∂θ_from = Base.RefValue{Float64}(0.0)
        ∂P∂V_from = Base.RefValue{Float64}(0.0)
        ∂Q∂V_from = Base.RefValue{Float64}(0.0)
        Vm_from = Vm[bus_from]
        for bus_to in data.neighbors[bus_from]
            if bus_to != bus_from
                col_to_vm = 2 * bus_to - 1
                col_to_va = 2 * bus_to
                bus_type = data.bus_type[bus_to, time_step]
                θ_from_to = θ[bus_from] - θ[bus_to]
                Vm_to = Vm[bus_to]
                Y_from_to = Yb[bus_from, bus_to]
                # 3 case if-else with Val(constant) is faster than 1 case with Val(bus_type)
                if bus_type == PSY.ACBusTypes.PQ
                    _set_entries_for_neighbor(Jv,
                        Y_from_to,
                        Vm_from,
                        Vm_to,
                        θ_from_to,
                        row_from_p,
                        row_from_q,
                        col_to_vm,
                        col_to_va,
                        ∂P∂θ_from,
                        ∂Q∂θ_from,
                        ∂P∂V_from,
                        ∂Q∂V_from,
                        Val(PSY.ACBusTypes.PQ))
                elseif bus_type == PSY.ACBusTypes.PV
                    _set_entries_for_neighbor(Jv,
                        Y_from_to,
                        Vm_from,
                        Vm_to,
                        θ_from_to,
                        row_from_p,
                        row_from_q,
                        col_to_vm,
                        col_to_va,
                        ∂P∂θ_from,
                        ∂Q∂θ_from,
                        ∂P∂V_from,
                        ∂Q∂V_from,
                        Val(PSY.ACBusTypes.PV))
                elseif bus_type == PSY.ACBusTypes.REF
                    _set_entries_for_neighbor(Jv,
                        Y_from_to,
                        Vm_from,
                        Vm_to,
                        θ_from_to,
                        row_from_p,
                        row_from_q,
                        col_to_vm,
                        col_to_va,
                        ∂P∂θ_from,
                        ∂Q∂θ_from,
                        ∂P∂V_from,
                        ∂Q∂V_from,
                        Val(PSY.ACBusTypes.REF))
                end
            end
        end
        col_from_vm = 2 * bus_from - 1
        col_from_va = 2 * bus_from
        # set entries in diagonal blocks
        if data.bus_type[bus_from, time_step] == PSY.ACBusTypes.PQ
            Jv[row_from_p, col_from_va] = ∂P∂θ_from[]
            Jv[row_from_q, col_from_va] = ∂Q∂θ_from[]
            ∂P∂V_from[] += 2 * real(Yb[bus_from, bus_from]) * Vm[bus_from]
            ∂Q∂V_from[] -= 2 * imag(Yb[bus_from, bus_from]) * Vm[bus_from]
            Jv[row_from_p, col_from_vm] = ∂P∂V_from[]
            Jv[row_from_q, col_from_vm] = ∂Q∂V_from[]
        elseif data.bus_type[bus_from, time_step] == PSY.ACBusTypes.PV
            Jv[row_from_q, col_from_vm] = -1.0
            Jv[row_from_p, col_from_va] = ∂P∂θ_from[]
            Jv[row_from_q, col_from_va] = ∂Q∂θ_from[]
        elseif data.bus_type[bus_from, time_step] == PSY.ACBusTypes.REF
            Jv[row_from_p, col_from_vm] = -1.0
            Jv[row_from_q, col_from_va] = -1.0
        end
    end
    return
end

"""
    calculate_loss_factors(data::ACPowerFlowData, Jv::SparseMatrixCSC{Float64, Int32}, time_step::Int)

Calculate and store the active power loss factors in the `loss_factors` matrix of the `ACPowerFlowData` structure for a given time step.

The loss factors are computed using the Jacobian matrix `Jv` and the vector `dSbus_dV_ref`, which contains the 
partial derivatives of slack power with respect to bus voltages. The function interprets changes in 
slack active power injection as indicative of changes in grid active power losses. 
KLU is used to factorize the sparse Jacobian matrix to solve for the loss factors.

# Arguments
- `data::ACPowerFlowData`: The data structure containing power flow information, including the `loss_factors` matrix.
- `Jv::SparseMatrixCSC{Float64, Int32}`: The sparse Jacobian matrix of the power flow system.
- `time_step::Int`: The time step index for which the loss factors are calculated.
"""
function calculate_loss_factors(
    data::ACPowerFlowData,
    Jv::SparseMatrixCSC{Float64, Int32},
    time_step::Int,
)
    bus_numbers = 1:first(size(data.bus_type))
    ref_mask = data.bus_type[:, time_step] .== (PSY.ACBusTypes.REF,)
    pvpq_mask = .!ref_mask
    ref = bus_numbers[ref_mask]
    pvpq = bus_numbers[pvpq_mask]
    pvpq_coords = Int32[]
    for i in pvpq
        push!(pvpq_coords, 2 * i - 1)  # 2x - 1
        push!(pvpq_coords, 2 * i)      # 2x
    end
    J_t = sparse(transpose(Jv[pvpq_coords, pvpq_coords]))
    dSbus_dV_ref = collect(Jv[2 .* ref .- 1, pvpq_coords])[:]
    fact = KLU.klu(J_t)
    lf = fact \ dSbus_dV_ref
    idx = 1:2:(2 * length(pvpq) - 1)  # only take the dPref_dP loss factors, ignore dPref_dQ
    data.loss_factors[pvpq_mask, time_step] .= lf[idx]
    data.loss_factors[ref_mask, time_step] .= 1.0
end
