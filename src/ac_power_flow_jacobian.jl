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
    # temp buffers for updating the Jacobian
    U::Vector{ComplexF64}
    V::Vector{ComplexF64}
    Ibus::Vector{ComplexF64}
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
    J.Jf!(J.Jv, J.data, time_step, J.U, J.V, J.Ibus)
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
    J.Jf!(J.Jv, J.data, time_step, J.U, J.V, J.Ibus)
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
    f = _update_jacobian_matrix_values!
    num_buses = size(get_bus_type(data), 1)
    return ACPowerFlowJacobian(
        data,
        f,
        Jv0,
        zeros(ComplexF64, num_buses),
        zeros(ComplexF64, num_buses),
        zeros(ComplexF64, num_buses),
    )
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

    numlines = length(get_branch_lookup(data))
    sizehint!(rows, 4 * numlines)
    sizehint!(columns, 4 * numlines)
    sizehint!(values, 4 * numlines)

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

"""Based on Algorithm 1 of \"Bus Admittance Matrix Revisited: Performance Challenges
on Modern Computers\"."""
function _update_jacobian_matrix_values!(
    Jv::SparseArrays.SparseMatrixCSC{Float64, Int32},
    data::ACPowerFlowData,
    time_step::Int64,
    U::Vector{ComplexF64},
    V::Vector{ComplexF64},
    Ibus::Vector{ComplexF64},
)
    Yb = data.power_network_matrix.data
    Vm = view(data.bus_magnitude, :, time_step)
    θ = view(data.bus_angles, :, time_step)

    copyto!(V, θ)
    V .*= im
    U .= exp.(V)
    copyto!(V, U)
    V .*= Vm
    Ibus .= zero(ComplexF64)
    # U = exp.(im .* θ)
    # V = Vm .* U
    # Ibus = zeros(ComplexF64, num_buses)

    num_buses = size(data.bus_type, 1)
    # WARNING: this reinterpret means indices from Jv.colptr or Jv.rowval are 2x too big!
    J_nonzeros = reinterpret(ComplexF64, SparseArrays.nonzeros(Jv))
    Yb_nonzeros = SparseArrays.nonzeros(Yb)
    # step 0: do equivalent of ∂S/∂θ .= Ybus, ∂S/∂V .= Ybus.
    J_colptr = Jv.colptr
    Yb_colptr = Yb.colptr
    Yb_rows = SparseArrays.rowvals(Yb)
    for col in 1:num_buses
        Yb_val_range = Yb.colptr[col]:(Yb_colptr[col + 1] - 1)
        J_Vm_start = div(J_colptr[2 * col - 1] + 1, 2)
        J_Vm_val_range = Yb_val_range .+ (J_Vm_start - Yb_val_range.start)
        Yb_col_view = @view Yb_nonzeros[Yb_val_range]
        J_Vm_view = @view J_nonzeros[J_Vm_val_range]
        copy!(J_Vm_view, Yb_col_view)
        J_Va_start = div(J_colptr[2 * col] + 1, 2)
        J_Va_val_range = Yb_val_range .+ (J_Va_start - Yb_val_range.start)
        J_Va_view = @view J_nonzeros[J_Va_val_range]
        copy!(J_Va_view, Yb_col_view)
    end
    # step 1: their "pass 1"
    for col in 1:num_buses
        Yb_val_range = Yb_colptr[col]:(Yb_colptr[col + 1] - 1)
        # Ibus[Yb_rows[Yb_val_range]] .+= Yb_nonzeros[Yb_val_range] .* V[col]
        for k in Yb_val_range
            Ibus[Yb_rows[k]] += Yb_nonzeros[k] * V[col]
        end

        J_Vm_start = div(J_colptr[2 * col - 1] + 1, 2)
        J_Vm_val_range = Yb_val_range .+ (J_Vm_start - Yb_val_range.start)
        J_nonzeros[J_Vm_val_range] .*= U[col]
        J_Va_start = div(J_colptr[2 * col] + 1, 2)
        J_Va_val_range = Yb_val_range .+ (J_Va_start - Yb_val_range.start)
        J_nonzeros[J_Va_val_range] .*= V[col]
    end
    # step 2: their "pass 2," with a little extra for the cols where input is power.
    for col in 1:num_buses
        bus_type = get_bus_type(data)[col, time_step]
        Yb_val_range = Yb_colptr[col]:(Yb_colptr[col + 1] - 1)
        bus_inds = @view Yb_rows[Yb_val_range]
        ind = searchsortedfirst(bus_inds, col)
        J_Vm_start = div(J_colptr[2 * col - 1] + 1, 2)
        J_Vm_val_range = Yb_val_range .+ (J_Vm_start - Yb_val_range.start)
        J_Vm_vals = @view J_nonzeros[J_Vm_val_range]
        if bus_type == PSY.ACBusTypes.PQ
            conj!(J_Vm_vals)
            # non-allocating version of J_Vm_vals .*= V[bus_inds]
            for i in Yb_val_range
                J_Vm_vals[i - Yb_val_range.start + 1] *= V[Yb_rows[i]]
            end
            J_nonzeros[(J_Vm_start - 1) + ind] += conj(Ibus[col]) * U[col]
        else
            J_Vm_vals .= zero(ComplexF64)
        end

        J_Va_start = div(J_colptr[2 * col] + 1, 2)
        J_Va_val_range = Yb_val_range .+ (J_Va_start - Yb_val_range.start)
        J_Va_vals = @view J_nonzeros[J_Va_val_range]
        if bus_type == PSY.ACBusTypes.PV || bus_type == PSY.ACBusTypes.PQ
            J_Va_vals .*= ComplexF64(-1)
            J_nonzeros[(J_Va_start - 1) + ind] += Ibus[col]
            conj!(J_Va_vals)
            # non-allocating version of J_Va_vals .*= im * V[bus_inds]
            for i in Yb_val_range
                J_Va_vals[i - Yb_val_range.start + 1] *= im * V[Yb_rows[i]]
            end
        else
            J_Va_vals .= zero(ComplexF64)
        end

        # correct for our power terms.
        if bus_type == PSY.ACBusTypes.REF
            J_nonzeros[(J_Vm_start - 1) + ind] = ComplexF64(-1)
            J_nonzeros[(J_Va_start - 1) + ind] = -im
        elseif bus_type == PSY.ACBusTypes.PV
            J_nonzeros[(J_Vm_start - 1) + ind] = -im
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
