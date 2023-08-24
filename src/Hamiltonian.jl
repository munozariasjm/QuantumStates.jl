@with_kw mutable struct Hamiltonian{T<:BasisState, F}
    basis::Vector{T}
    operator::Expr
    parameters::ParameterList
    states::Vector{State{T}}           = states_from_basis(basis)
    operators::F                       = unpack_operator(operator, basis)
    matrix::Matrix{ComplexF64}         = sum(parameters.param_dict[operator.param] .* operator.matrix for operator ∈ operators)
    tdms::Array{ComplexF64, 3}         = zeros(ComplexF64, length(states), length(states), 3)
    tdms_m::Array{ComplexF64, 3}       = zeros(ComplexF64, length(states), length(states), 3)
    basis_tdms::Array{ComplexF64, 3}   = zeros(ComplexF64, length(basis), length(basis), 3)
    basis_tdms_m::Array{ComplexF64, 3} = zeros(ComplexF64, length(basis), length(basis), 3)
end
export Hamiltonian

"""
    
"""
@with_kw mutable struct CombinedHamiltonian{T1<:BasisState, T2<:BasisState, F, F1, F2}
    H1::Hamiltonian{T1, F1}
    H2::Hamiltonian{T2, F2}
    P::Matrix{ComplexF64}           = change_of_basis_matrix(H1.basis, H2.basis)
    basis::Vector{T1}               = H1.basis
    operator::Expr                  = Expr(:nothing)
    parameters::ParameterList       = ParameterList(Dict{Symbol, Float64}())
    states::Vector{State{T1}}       = [H1.states; convert_basis(H2.states, basis)]
    operators::F                    = unpack_operator(operator, basis)
    matrix::Matrix{ComplexF64}      = create_block_diagonal_matrix(H1.matrix, H2.matrix)
    tdms::Array{ComplexF64, 3}      = zeros(ComplexF64, length(states), length(states), 3)
    basis_tdms::Array{ComplexF64}   = zeros(ComplexF64, length(basis), length(basis), 3)
end
export CombinedHamiltonian

function change_of_basis_matrix(basis, basis′)
    P = zeros(ComplexF64, length(basis), length(basis′))
    for (i, state) in enumerate(basis)
        for (j, state′) in enumerate(basis′)
            P[i,j] = overlap(state, state′)
        end
    end
    return P
end

function convert_basis(H::Hamiltonian, basis′)

    states′ = State{typeof(basis′[1])}[]
    basis = H.states[1].basis
    states = H.states

    P = change_of_basis_matrix(basis, basis′)

    for state ∈ states
        state′ = State(state.E, basis′, P' * state.coeffs, state.idx)
        push!(states′, state′)
    end

    matrix′ = P' * H.matrix * P

    operators′ = Operator[]
    for operator ∈ H.operators
        operator_matrix′ = P' * operator.matrix * P
        push!(operators′, Operator(operator.param, operator.operator, operator_matrix′))
    end

    return Hamiltonian(
        basis=basis′,
        operator=H.operator,
        parameters=H.parameters,
        states=states′,
        operators=operators′,
        matrix=matrix′,
        tdms=H.tdms,
        tdms_m=H.tdms_m,
        basis_tdms=zeros(ComplexF64, length(basis′), length(basis′), 3)
    )
end
export convert_basis

function extend_basis!(H::Hamiltonian, basis′)
    extended_basis = unique([H.basis; basis′])
    size_increase = length(extended_basis) - length(H.basis)
    H.basis = extended_basis
    for i ∈ eachindex(H.states)
        H.states[i].basis = extended_basis
        H.states[i].coeffs = vcat(H.states[i].coeffs, zeros(size_increase))
    end
    return nothing
end
export extend_basis!

# """
#     Assumes that `H₁` and `H₂` have the same operator and set of parameters. 
# """
# function combine(H₁::Hamiltonian, H₂::Hamiltonian)
#     combined_basis = [H₁.basis; H₂.basis]
#     combined_states = [H₁.states; H₂.states]
#     Hamiltonian(basis=combined_basis, operator=H₁.operator, parameters=H₁.parameters)
# end

subspace(H::Hamiltonian, QN_bounds) = Hamiltonian(basis=H.basis, operator=H.operator, parameters=H.parameters, states=subspace(H.states, QN_bounds)[2])
export subspace
# subspace(H::Hamiltonian, basis_idxs) = Hamiltonian(H.basis[basis_idxs], H.operator, H.parameters)

function update_basis_tdms!(H::Hamiltonian)
    for i ∈ eachindex(H.basis)
        for j ∈ (i+1):length(H.basis)
            bstate = H.basis[i]
            bstate′ = H.basis[j]
            for p ∈ -1:1
                H.basis_tdms[i,j,p+2] = TDM(bstate, bstate′, p)
                H.basis_tdms[j,i,p+2] = conj(H.basis_tdms[i,j,p+2])
            end
        end
    end
    return nothing
end
export update_basis_tdms!

function update_basis_tdms_m!(H::Hamiltonian)
    for (i, bstate) ∈ enumerate(H.basis)
        for (j, bstate′) ∈ enumerate(H.basis)
            for p ∈ -1:1
                H.basis_tdms_m[i,j,p+2] = TDM_magnetic(bstate, bstate′, p)
                # H.basis_tdms_m[j,i,p+2] = conj(H.basis_tdms_m[i,j,p+2])
            end
        end
    end
    return nothing
end
export update_basis_tdms_m!

# Updated 7/15/2023
function update_tdms!(H::Hamiltonian, idxs=eachindex(H.states))
    for i ∈ idxs, j ∈ idxs
        if j >= i
            state = H.states[i]
            state′ = H.states[j]
            for p ∈ -1:1
                H.tdms[i,j,p+2] = H.tdms[j,i,p+2] = 0.0
                for m ∈ eachindex(H.basis), n ∈ eachindex(H.basis)
                    H.tdms[i,j,p+2] += conj(state.coeffs[m]) * state′.coeffs[n] * H.basis_tdms[m,n,p+2]
                end
                H.tdms[j,i,p+2] = conj(H.tdms[i,j,p+2])
            end
        end
    end
    return nothing
end
export update_tdms!

# function update_tdms!(H::Hamiltonian, idxs=eachindex(H.states))
#     for i ∈ idxs, j ∈ idxs
#         if j >= i
#             state = H.states[i]
#             state′ = H.states[j]
#             for p ∈ -1:1
#                 H.tdms[i,j,p+2] = H.tdms[j,i,p+2] = 0.0
#                 for m ∈ eachindex(H.basis), n ∈ eachindex(H.basis)
#                     H.tdms[i,j,p+2] += conj(state.coeffs[m]) * state′.coeffs[n] * H.basis_tdms[m,n,p+2]
#                 end
#                 H.tdms[j,i,p+2] += conj(H.tdms[i,j,p+2])
#             end
#         end
#     end
#     return nothing
# end
# export update_tdms!

function update_tdms_m!(H::Hamiltonian, idxs=eachindex(H.states))
    for i ∈ idxs, j ∈ idxs
        if j >= i
            state = H.states[i]
            state′ = H.states[j]
            for p ∈ -1:1
                H.tdms[i,j,p+2] = H.tdms[j,i,p+2] = 0.0
                for m ∈ eachindex(H.basis), n ∈ eachindex(H.basis)
                    H.tdms[i,j,p+2] += conj(state.coeffs[m]) * state′.coeffs[n] * H.basis_tdms_m[m,n,p+2]
                end
                H.tdms[j,i,p+2] += conj(H.tdms[i,j,p+2])
            end
        end
    end
    return nothing
end
export update_tdms_m

function add_to_H(H::Hamiltonian, param::Symbol, f::Function)
    operator = Expr(:call, :+, H.operator, f)
    matrix = matrix_from_operator(H.basis, f)
    param_dict = Dict(H.parameters.param_dict..., param => 0.0) # any new term added has its parameter value set to zero as default
    parameters = ParameterList(param_dict)
    operators = (; H.operators..., param => Operator(param, f, matrix))
    return Hamiltonian(H.basis, operator, parameters, H.states, operators, H.matrix, H.tdms, H.tdms_m, H.basis_tdms, H.basis_tdms_m)
end
export add_to_H

function add_to_H(H::CombinedHamiltonian, param::Symbol, f::Function)
    operator = Expr(:call, :+, H.operator, f)
    matrix = matrix_from_operator(H.basis, f)
    param_dict = Dict(H.parameters.param_dict..., param => 0.0) # any new term added has its parameter value set to zero as default
    parameters = ParameterList(param_dict)
    operators = (; H.operators..., param => Operator(param, f, matrix))

    return CombinedHamiltonian(
        H1=H.H1,
        H2=H.H2,
        basis=H.basis,
        operator=operator,
        parameters=parameters,
        states=H.states,
        operators=operators,
        matrix=H.matrix,
        tdms=H.tdms,
        basis_tdms=H.basis_tdms
    )
end
export add_to_H

# function +(H, O::Operator)
#     H.operator = Expr(:call, :+, H.operator, O)
#     H.parameters.param_dict[O.param] = 0.0 # any new term added has its parameter value set to zero as default
#     H.operators = (; H.operators..., O.operator)
#     print(H.parameters)
#     return H
# end

# function +(H, expr::Expr)
#     param = expr.args[2]
#     operator = eval(expr.args[3])
#     matrix = matrix_from_operator(H.basis, operator)
#     O = Operator(param, operator, matrix)
#     return H + O
# end
    
# function +(H1::Hamiltonian, H2::Hamiltonian)
#     Hamiltonian(H1.basis, )


@generated function evaluate_(H::Hamiltonian, _::Val{N}) where N
    quote
        Base.Cartesian.@nexprs $N i -> evaluate_operator!(H.basis, H.operators[i])
    end
end
export evaluate_

function add_operator_matrix!(H, operator::Operator)
    operator_val = getproperty(H.parameters, operator.param)
    for i ∈ eachindex(H.matrix)
        H.matrix[i] += operator_val * operator.matrix[i]
    end
    return nothing
end
export add_operator_matrix!

@generated function sum_operator_matrices!(H, _::Val{N}) where N
    quote
        Base.Cartesian.@nexprs $N i -> add_operator_matrix!(H, H.operators[i])
    end
end

"""
    full_evaluate!(H::Hamiltonian)

    Additionally reevaluates the matrix elements of all operators in the Hamiltonian. 
"""
function full_evaluate!(H::Hamiltonian)
    evaluate_(H, Val(length(H.operators)))
    H.matrix = sum(getproperty(H.parameters, operator.param) .* operator.matrix for operator ∈ H.operators)
    return nothing
end
export full_evaluate!

function evaluate!(H::Hamiltonian)
    H.matrix .= 0.0
    sum_operator_matrices!(H, Val(length(H.operators)))
    return nothing
end

function evaluate!(H::CombinedHamiltonian)
    H.matrix .= 0.0
    if length(H.operators) > 0
        sum_operator_matrices!(H, Val(length(H.operators)))
    end
    evaluate!(H.H1)
    evaluate!(H.H2)
    n = length(H.H1.states)
    m = length(H.H2.states)
    H.matrix[1:n, 1:n] .+= H.H1.matrix
    size(H.P') |> display
    size(H.H2.matrix) |> display
    size(H.P) |> display
    # size(H.P' * H.H2.matrix * H.P') |> display
    print(m)
    # H.matrix[(n+1):(n+m), (n+1):(n+m)] .+= H.P' * H.H2.matrix * H.P
    return nothing
end
export evaluate!
        
function solve!(H)
    es, vs = eigen(H.matrix)
    for i ∈ eachindex(H.states)
        H.states[i].E = real(es[i])
        H.states[i].coeffs = vs[:,i]
    end
    return nothing
end
export solve!

# function solve!(H::CombinedHamiltonian)
#     solve!(H.H1)
#     solve!(H.H2)
#     H.states = [H.H1.states; convert_basis(H.H2.states, H.basis)]
#     return nothing
# end
# export solve!

@with_kw mutable struct Hamiltonian_Old{T1, T2<:BasisState}
    H_operator::T1
    basis::Vector{T2}
    states::Vector{State{T2}} = states_from_basis(basis)
    M::Array{Float64, 2}      = zeros(Float64, length(basis), length(basis))
    M_tdms::Array{Float64, 2} = zeros(Float64, length(basis), length(basis))
#     M_basis::Array{Float64, 2}  = zeros(Float64, length(basis), length(basis))
#     M_states::Array{Float64, 2} = reduce(vcat, [state.coeffs for state in states]')
end
export Hamiltonian_Old

function energies(H::Hamiltonian_Old)
    return [state.E for state ∈ H.states]
end
export energies

# import LinearAlgebra: mul!
# function mul!(H_new, H::Hamiltonian, M)
# #     H_new = zeros(Float64, size(M,1), size(M,2))
#     for n ∈ indices(M,2), m ∈ eachindex(H.states)
#         coeffs = H.states[m].coeffs
#         Cmn = zero(eltype(M))
#         for k ∈ eachindex(coeffs)
#             Cmn += coeffs[k] * M[k,n]
#         end
#         H_new[m,n] = Cmn
#     end
#     return nothing
# end
# export mul!

# function *(H::Hamiltonian, M)
#     H_new = zeros(Float64, size(M,1), size(M,2))
#     for n ∈ indices(M,2), m ∈ eachindex(H.states)
#         coeffs = H.states[m].coeffs
#         Cmn = zero(eltype(M))
#         for k ∈ eachindex(coeffs)
#             Cmn += coeffs[k] * M[k,n]
#         end
#         H_new[m,n] = Cmn
#     end
#     return H_new
# end
# export +, -, *, ^
        
function update_matrix!(H::Hamiltonian_Old, p, args::Vararg{Any, T}) where {T}
    basis = H.basis
    for (i, state) ∈ enumerate(basis)
        for (j, state′) ∈ enumerate(basis)
            H.M[i,j] = H.H_operator(state, state′, i, j, p, args...)
        end
    end
    return nothing
end
export update_matrix!

function solve!(H::Hamiltonian_Old)
    es, vs = eigen(H.M)
    for i in eachindex(H.states)
        H.states[i].E = real(es[i])
        H.states[i].coeffs = real(vs[:,i])
    end
    return nothing
end
export solve!

# function collapse!(H::Hamiltonian, collapsed_QNs)
#     """
#     Reduces a Hamiltonian by aggregating states that have the same quantum numbers,
#     other than those included in `collapsed_QNs`. For each one of the aggregated 
#     states, the corresponding entry in the Hamiltonian is the sum of all the entries
#     for the original states that "collapse" into this state.
#     """
#     state_type = typeof(H.basis[1])
#     QNs = fieldnames(state_type)[1:end-1]
    
#     new_QNs = setdiff(QNs, collapsed_QNs)
#     new_basis = BasisState[]
#     basis_mapping = zeros(Int64, length(H.basis))
#     prev_QN_vals = zeros(Rational, length(new_QNs))
#     prev_QN_vals .= Inf
#     idx = 0
#     multi
    
#     for (i, state) in enumerate(H.basis)
#         multiplicity = 0
#         added_state = false
#         for (j, QN) in enumerate(new_QNs)
#             QN_val = getfield(state, QN)
#             if QN_val != prev_QN_vals[j] && !added_state
#                 push!(new_basis, state)
#                 added_state = true
#                 idx += 1
#             end
#             prev_QN_vals[j] = QN_val
#         end
#         basis_mapping[i] = idx
#     end
    
#     new_H = zeros(Complex, length(new_basis), length(new_basis))
#     for (i, state) in enumerate(H.basis)
#         for (j, state) in enumerate(H.basis)
#             idx1 = basis_mapping[i]
#             idx2 = basis_mapping[j]
#             new_H[idx1,idx2] += H.M[i,j]
#         end
#     end
    
#     H.M = new_H
#     H.basis = new_basis
#     return nothing
# end
# export collapse!

# function solve(H::Hamiltonian; kwargs...)
#     basis = H.basis
#     H_basis = zeros(ComplexF64, length(basis), length(basis))

#     args = [kwarg.second for kwarg in kwargs]
#     if !isempty(args)
#         H_operator = H.H_operator(args...)
#     else
#         H_operator = H.H_operator
#     end
#     for (i, state) in enumerate(basis)
#         for (j, state′) in enumerate(basis)
#             H_basis[i,j] = H_operator(state, state′)
# #             println(H_basis[i,j])
#         end
#     end
#     H.M = H_basis #states * H_basis * states'
# #     eigvals, eigvecs = eigen(H.M)
#     return eigen(H.M)
# end
# export solve

### Parameter scanning
# @with_kw mutable struct ParameterScan2{N, T}
#     scan_params::NTuple{N, Symbol}
#     state_dict::Dict{NTuple{N, Float64}, State{T}}
# end
# export ParameterScan2

@with_kw mutable struct ParameterScan1{T<:BasisState}
    state_dict::Dict{Float64, Vector{State{T}}}
end

# import Base: getindex
# function getindex(parameter_scan::ParameterScan{T}, params) where {T}
#     H_dict = parameter_scan.H_dict
#     return H_dict[(params[scan_param] for scan_param ∈ parameter_scan.scan_params)...]
# end

function scan_single_parameter(H::Hamiltonian, param::Symbol, scan_range)
    
    d = Dict{Float64, Vector{State{typeof(H.basis[1])}}}()
    parameter_scan = ParameterScan1(d)
    overlaps = zeros(Float64, length(H.states), length(H.states))

    prev_states = deepcopy(H.states)
    for (i, scan_value) ∈ enumerate(scan_range)
        
        H.parameters.param_dict[param] = scan_value
        evaluate!(H)
        solve!(H)
        
        # Ensure that state indices are correctly tracked, i.e., crossing states are not swapped
        if i != 1
            calculate_state_overlaps!(H.states, prev_states, overlaps)
            new_idxs = tracked_idxs(overlaps)
            H.states = H.states[new_idxs]

            for j in eachindex(H.states)
                # H.states[j].idx = j
                prev_states[j].coeffs = H.states[j].coeffs
            end
        end
        parameter_scan.state_dict[scan_value] = deepcopy(H.states)
        
    end
    return parameter_scan
end
export scan_single_parameter
export scan_parameters

"""
    scan_parameters()

    scan_iterator = 
"""
function scan_parameters(H::Hamiltonian, scan_values::T, iterator::F1, H_func!::F2, output_func::F3; n_threads=Threads.nthreads()) where {T,F1,F2,F3}
    scan_iterator = iterator(scan_values...)
    
    n_values = length(scan_iterator)
    # saved_values = zeros(n_values, length(H.states))

    output_type = typeof(output_func(H))

    n_params = length(first(scan_iterator))
    # print(n_params)

    prog_bar = Progress(n_values)

    # Multi-threading settings
    batch_size = cld(n_values, n_threads)
    remainder = n_values - batch_size * n_threads
    partitions = Iterators.partition(scan_iterator, batch_size)
    # println(length.(partitions))

    # full_dict = Dict{NTuple{n_params, Float64}, output_type}()

    # print(length(partitions))
    tasks = Vector{Task}(undef, n_threads)

    H_copy = deepcopy(H)
    H_func!(H_copy, scan_values[1])

    @sync for (i, scan_values) ∈ enumerate(partitions)

        _H = deepcopy(H_copy)
        tracked_idxs = collect(1:length(_H.states))
        overlaps = zeros(Float64, length(_H.states), length(_H.states))
        prev_states = deepcopy(_H.states)
        dict = Dict{NTuple{n_params, Float64}, output_type}()
        dict_state_idxs = Dict{NTuple{n_params, Float64}, Vector{Int64}}()
        state_idxs = collect(1:length(_H.states))

        @async tasks[i] = Threads.@spawn begin
        # _batch_size = i <= remainder ? (batch_size + 1) : batch_size - 1
        # batch_start_idx = 1 + ((i <= remainder) ? i : remainder) + batch_size * (i-1)
        # for j ∈ batch_start_idx:(batch_start_idx + _batch_size)
            for scan_value ∈ scan_values
                H_func!(_H, scan_value)
                # if j != 1
                calculate_state_overlaps!(_H.states, prev_states, overlaps)
                tracked_idxs!(overlaps, _H.states, prev_states, tracked_idxs)
                _H.states .= _H.states[tracked_idxs]
    
                for j ∈ eachindex(_H.states)
                    prev_states[j].E = _H.states[j].E
                    prev_states[j].coeffs .= _H.states[j].coeffs
                end
                # end
                # saved_values[i,:] = output_func(H)
                dict[scan_value] = output_func(_H)
                # state_idxs .= state_idxs[tracked_idxs]
                dict_state_idxs[scan_value] = deepcopy(state_idxs[tracked_idxs])
                next!(prog_bar)
            end
            dict, dict_state_idxs
        end
        # full_dict = merge(full_dict, fetch(task))
    end

    # Returned results are a vector of tuples holding dictionaries (with energies) and vectors of state indices
    scan_data = fetch.(tasks)
    # print([x[1] for x ∈ scan_data][2])
    full_dict = merge([x[1] for x ∈ scan_data]...)
    full_dict_tracked_idxs =  merge([x[2] for x ∈ scan_data]...)

    # for (i, scan_values) ∈ enumerate(scan_iterator)
        
    #     _H = deepcopy(H)
    #     H_func!(H, scan_values, i)
        
    #     # Ensure that state indices are correctly tracked, i.e., crossing states are not swapped
    #     if i != 1
    #         calculate_state_overlaps!(H.states, prev_states, overlaps)
    #         tracked_idxs!(overlaps, tracked_idxs)
    #         H.states = H.states[tracked_idxs]

    #         for j ∈ eachindex(H.states)
    #             prev_states[j].coeffs .= H.states[j].coeffs
    #         end
    #     end
    #     # saved_values[i,:] = output_func(H)
    #     dict[scan_values] = output_func(H)
    #     next!(prog_bar)
    # end
    return sort(full_dict), sort(full_dict_tracked_idxs)
    
end
export scan_parameters

# function scan_parameters(H::Hamiltonian, scan_ranges) where {T}
    
#     scan_params = keys(scan_ranges)
    
#     H_dict = Dict{NTuple{length(scan_params), Float64}, Vector{State}}()
#     parameter_scan = ParameterScan(scan_params=scan_params, H_dict=H_dict)
    
#     scan_params = keys(scan_ranges)
#     scan_ranges_product = Iterators.product(values(scan_ranges)...)
    
#     prev_states = deepcopy(H.states)
# #     idxs_prev = collect(1:length(H.basis))
#     for (i, scan_values) ∈ enumerate(scan_ranges_product)
        
#         H.scan_param
#         update_matrix!(H_scan, scan_values...)
#         solve!(H_scan)
        
#         # Ensure that state indices are correctly tracked, i.e., crossing states are not swapped
#         if i != 1
#             overlaps = calculate_state_overlaps(H_scan.states, H_prev.states)
#             new_idxs = tracked_idxs(overlaps)
            
#             H_scan.states = H_scan.states[new_idxs]
#             for j in eachindex(H_scan.states)
#                 H_scan.states[j].idx = j
#             end
            
#             H_prev = H_scan
# #             idxs_prev = idxs_prev[new_idxs]
#         end
        
#         parameter_scan.H_dict[scan_values] = H_scan
        
#     end
#     return parameter_scan
# end
# export scan_parameters

import Plots: plot
function plot(parameter_scan::ParameterScan1, plot_ranges, plot_parameter; subspace_idxs=nothing)
    
    H_dict = parameter_scan.H_dict
    
    scan_ranges = collect(keys(H_dict))
    H_any = H_dict[scan_ranges[1]]
    
    xs = Float64[]
    
    if isnothing(subspace_idxs)
        ys = Array{Float64, 2}(undef, 0, length(H_any.basis))
    else
        ys = Array{Float64, 2}(undef, 0, length(subspace_idxs))
    end
    
    plot_ranges_product = Iterators.product(values(plot_ranges)...)
    
    plot_parameter_idx = -1
    for (i, scan_param) ∈ enumerate(parameter_scan.scan_params)
         if scan_param == plot_parameter
            plot_parameter_idx = i
        end
    end
    
    for plot_values ∈ plot_ranges_product
        x = plot_values[plot_parameter_idx]
        push!(xs, x)
        
        H = H_dict[plot_values]
        es = energies(H)
        
        if isnothing(subspace_idxs)
            ys = [ys; es']
        else
            ys = [ys; es[subspace_idxs]']
        end
    end
        
    plot(xs, ys, legend=nothing)
end
export plot

"""
Implement functionality to save predefined Hamiltonians.
"""

using Serialization

function save_to_file(H::Hamiltonian, name::String, filepath::String)
    serialize(filepath * name * ".jl", H)
end
export save_to_file

function load_from_file(name::String, filepath::String)
    deserialize(filepath * name * ".jl")
end
export load_from_file
                                                        
# using JLD2

# function save_to_file(H::Hamiltonian, name::String, filepath::String)
#     save(filepath * name * ".jld2", name, H)
# end
# export save_to_file

# function load_from_file(name::String, filepath::String)
#     load(filepath * name * ".jld2", name)
# end
# export load_from_file

### Only finds the states where the QNs are good quantum numbers
# function subspace(states::Vector{State}, QN_bounds)
#     subspace = State[]
#     QNs = keys(QN_bounds)
#     add_to_subspace = ones(Bool, length(states))
#     for QN ∈ QNs
#         for (i, state) ∈ enumerate(states)
#             is_QN_good, QN_vals = is_good_quantum_number(state, QN)
#             if !is_QN_good
#                 error("Quantum number is not well-defined.")
#             elseif !(QN_vals[1] ∈ QN_bounds[QN])
#                 add_to_subspace[i] = false
#             end
#         end
#     end
#     for i ∈ eachindex(states)
#         if add_to_subspace[i]
#             push!(subspace, states[i])
#         end
#     end
#     return subspace
# end
# export subspace