# import WignerSymbols_Simple: wigner3j, wigner6j
import WignerSymbols
import LinearAlgebra: norm

using CompositeStructs
using DataFrames
using PartialWaveFunctions

using HalfIntegers
using Parameters
using NamedTupleTools
using LinearAlgebra
using LoopVectorization
using NamedTupleTools

import ProgressMeter: Progress, next!

abstract type BasisState end
export BasisState

Base.@kwdef struct QuantumState <: BasisState
    E::Float64 = 0.0
end

mutable struct State{T<:BasisState}
    E::Float64
    basis::Vector{T}
    coeffs::Vector{ComplexF64}
    idx::Int64
end
export State

⊕(f1::T1, f2::T2) where {T1,T2} = (bs, bs′) -> f1(bs.basis_state1, bs′.basis_state1) + f2(bs.basis_state2, bs′.basis_state2)
export ⊕

⊗(f1::T1, f2::T2) where {T1,T2} = (bs, bs′) -> f1(bs.basis_state1, bs′.basis_state1) * f2(bs.basis_state2, bs′.basis_state2)
export ⊗

function expectation(state::State, s::Symbol)
    exp = zero(ComplexF64)
    for i ∈ eachindex(state.basis)
        exp += getfield(state.basis[i], s) * state.coeffs[i] * conj(state.coeffs[i])
    end
    return exp
end
export expectation 

energy(s::State) = s.E
export energy

"""
Returns a vector of states based on a basis.
"""
function states_from_basis(basis::Vector{<:BasisState})
    states = [State(0.0, basis, zeros(ComplexF64, length(basis)), i) for i in 1:length(basis)]
    for i in eachindex(basis)
        states[i].coeffs[i] = 1.0
    end
    return states
end
export states_from_basis

function states_to_basis(states::Vector{State{T}}) where T
    basis = states[1].basis
    basis_states_include = zeros(Bool, length(basis))
    for i ∈ eachindex(states)
        state = states[i]
        normed_state_coeffs = norm.(state.coeffs)
        nonzero_idxs = findall(x -> x > 0, normed_state_coeffs)
        basis_states_include[nonzero_idxs] .|= true
    end
    return findall(>(false), basis_states_include), basis[basis_states_include]
end
export states_to_basis

import LinearAlgebra: dot
dot(state::State, state′::State) = state.coeffs ⋅ state′.coeffs
export ⋅

norm(state::State) = norm(state.coeffs)
export norm

function wigner3j_fromcg(j1, m1, j2, m2, j3, m3)
    return (-1)^(j1-j2-m3) * CG(j1, m1, j2, m2, j3, -m3) / sqrt(2j3 + 1)
end
export wigner3j_fromcg

function wigner6j_fromcg(j1, j2, j3, j4, j5, j6)
    sum(
        (-1)^(j1 + j2 + j3 + j4 + j5 + j6 - m1 - m2 - m3 - m4 - m5 - m6)
        * wigner3j_fromcg(j1, -m1, j2, -m2, j3, -m3)
        * wigner3j_fromcg(j1, m1, j5, -m5, j6, m6)
        * wigner3j_fromcg(j4, m4, j2, m2, j6, -m6)
        * wigner3j_fromcg(j4, -m4, j5, m5, j3, m3)
        for m1 ∈ -j1:j1,
            m2 ∈ -j2:j2,
            m3 ∈ -j3:j3,
            m4 ∈ -j4:j4,
            m5 ∈ -j5:j5,
            m6 ∈ -j6:j6
        if (m1 + m2 + m3 == 0) & (m1 - m5 + m6 == 0) & (m4 + m2 - m6 == 0) & (-m4 + m5 + m3 == 0)
    )
end
export wigner6j_fromcg

function wigner3j_(j1, j2, j3, m1, m2, m3)
    try 
        WignerSymbols.wigner3j(Float64, j1, j2, j3, m1, m2, m3) 
    catch 
        0.0
    end
end

function wigner6j_(j1, j2, j3, m1, m2, m3)
    try 
        WignerSymbols.wigner6j(Float64, j1, j2, j3, m1, m2, m3) 
    catch 
        0.0
    end
end
export wigner6j_

function wigner9j(j1, j2, j3, j4, j5, j6, j7, j8, j9)::Float64
    val = 0.0
    kmin = max(abs(j1 - j9), abs(j4 - j8), abs(j2 - j6))
    kmax = min(abs(j1 + j9), abs(j4 + j8), abs(j2 + j6))
    if kmax >= kmin
        val += sum(
            (-1)^(2k) * (2k + 1) * 
            wigner6j(j1, j4, j7, j8, j9, k) * 
            wigner6j(j2, j5, j8, j4, k, j6) *
            wigner6j(j3, j6, j9, k, j1, j2) for k in kmin:kmax)
    end 
    return val
end

function order(basis::Vector{<:BasisState}, ordering_QN)
    all_QN = [getfield(state, ordering_QN) for state in basis]
    idxs_sorted = sortperm(all_QN)
    return basis[idxs_sorted]
end
export order

function is_good_quantum_number(state::State, QN::Symbol)
    vals = Rational[]
    for (i, coeff) in enumerate(state.coeffs)
        val = getfield(state.basis[i], QN)
        if norm(coeff) > 1e-5
            push!(vals, val)
        end
    end
    unique_vals = unique(vals)
    return (length(unique_vals) == 1, unique_vals)
end

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

function subspace(states::Vector{State{T}}, QN_bounds, threshold=0.01) where {T}
    subspace = State{T}[]
    subspace_idxs = Int64[]
    QNs = keys(QN_bounds)
    add_to_subspace = ones(Bool, length(states))
    for QN ∈ QNs
        for (i, state) ∈ enumerate(states)
            for (j, coeff) ∈ enumerate(state.coeffs)
                if getfield(state.basis[j], QN) ∉ QN_bounds[QN]
                    if norm(coeff)^2 > threshold
                        add_to_subspace[i] = false
                    end
                end
            end
        end
    end
    for i ∈ eachindex(states)
        if add_to_subspace[i]
            push!(subspace, states[i])
            push!(subspace_idxs, i)
        end
    end
    return (subspace_idxs, subspace)
end

function subspace(basis::Vector{<:BasisState}, QN_bounds)
    subspace_basis = typeof(basis[1])[]
    QNs = keys(QN_bounds)
    add_to_subspace = ones(Bool, length(basis))
    for QN ∈ QNs
        for (i, b_state) ∈ enumerate(basis)
            if getfield(b_state, QN) ∉ QN_bounds[QN]
                add_to_subspace[i] = false
            end
        end
    end
    for i ∈ eachindex(basis)
        if add_to_subspace[i]
            push!(subspace_basis, basis[i])
        end
    end
    return subspace_basis
end
export subspace

function extend_operator(operator::T, state::State, state′::State, args...) where {T}
    val = zero(ComplexF64)
    for (i, basis_state) in enumerate(state.basis)
        for (j, basis_state′) in enumerate(state′.basis)
            val += state.coeffs[i] * state′.coeffs[j] * operator(basis_state, basis_state′, args...)
            # coeff1 = state.coeffs[i]
            # coeff2 = state′.coeffs[i]
            # if (norm(coeff1)^2 > 1e-5) && (norm(coeff1)^2 > 1e-5)
                # val += state.coeffs[i] * state′.coeffs[j] * operator(basis_state, basis_state′, args...)
            # end
        end
    end
    return val
end
export extend_operator

function TDM(state::State, state′::State, args...)
    tdm = zero(ComplexF64)
    for (i, basis_state) in enumerate(state.basis)
        for (j, basis_state′) in enumerate(state′.basis)
            tdm += conj(state.coeffs[i]) * state′.coeffs[j] * TDM(basis_state, basis_state′, args...)
            # coeff1 = state.coeffs[i]
            # coeff2 = state′.coeffs[i]
            # if (norm(coeff1)^2 > 1e-5) && (norm(coeff2)^2 > 1e-5)
                # tdm += state.coeffs[i] * state′.coeffs[j] * TDM(basis_state, basis_state′, args...)
            # end
        end
    end
    return tdm
end
export TDM

function convert_basis(states::Vector{State{T}}, basis′) where {T}

    new_states = State{typeof(basis′[1])}[]
    basis = states[1].basis

    # Calculate change-of-basis matrix
    P = zeros(Float64, length(basis′), length(basis))
    for (i, state) in enumerate(basis′)
        for (j, state′) in enumerate(basis)
            P[i,j] = overlap(state, state′)
        end
    end

    for state in states
        new_state = State(state.E, basis′, P * state.coeffs, state.idx)
        push!(new_states, new_state)
    end
    
#     # Check that the basis provided included all necessary states to make the conversion
#     for state in new_states
#         if norm(state.coeffs) 
    
    return new_states
end
export convert_basis

# function convert_basis(H::Hamiltonian, basis′)

#     P = deepcopy(H.matrix)
#     P .= 0.0

#     basis = H.basis
#     for (i, state) in enumerate(basis)
#         for (j, state′) in enumerate(basis′)
#             P[i,j] = overlap(state, state′)
#         end
#     end
#     H.M .= P' * H.M * P
#     H.basis = basis′

#     return Hamiltonian()
# end

# function minimize_basis(states::Vector{State})
#     basis_idxs = zeros(Bool, length(states[1].basis))
#     for state in states
#         for (i, coeff) in enumerate(state.coeffs)
#             if norm(coeff)^2 > 1e-5
#                 basis_idxs[i] = true
#             end
#         end
#     end
#     return

function states_table(states::Vector{<:State}; threshold=1e-3, dominant_state_only=false)
    QNs = fieldnames(typeof(states[1].basis[1]))[2:end-1]
    QN_columns = NamedTuple(QN => Rational{Int64}[] for QN in QNs)
    all_columns = merge((State = Integer[], c = Float64[],), QN_columns, (E = Float64[],))
    df = DataFrame(all_columns)

    vals = zeros(Float64, length(QNs)+3)
    for (i, state) ∈ enumerate(states)
        basis_states, coeffs = contributing_basis_states(state, threshold)
        vals[1] = i
        vals[end] = state.E
        
        cs = norm.(coeffs).^2
        if dominant_state_only
            max_c, dominant_state_idx = findmax(cs)
            basis_state = basis_states[dominant_state_idx]
            vals[2] = max_c
            for (i, QN) ∈ enumerate(QNs)
                vals[i+2] = getfield(basis_state, QN)
            end
            push!(df, vals)
        else
            for (basis_state, coeff) ∈ zip(basis_states, coeffs)
                c = norm(coeff)^2
                vals[2] = c
                for (i, QN) ∈ enumerate(QNs)
                    vals[i+2] = getfield(basis_state, QN)
                end
                push!(df, vals)
            end
        end
    end
    return df
end

states_table(states::Vector{<:State}, relabelling_states::Vector{<:State}; threshold=1e-3, dominant_state_only=false) =
    states_table(relabelling_states, threshold=threshold, dominant_state_only=dominant_state_only)
export states_table
        
# Extend the base "+" function such that we can simply add matrix elements before calling them
import Base: +, -, *, /, ^

+(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) + g(args...; kwargs...)
-(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) - g(args...; kwargs...)
*(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) * g(args...; kwargs...)
*(c::Number, f::Function) =
    (args...; kwargs...) -> c * f(args...; kwargs...)
/(f::Function, c::Number) =
    (args...; kwargs...) -> f(args...; kwargs...) / c
^(f::Function, c::Real) = 
    (args...; kwargs...) -> f(args...; kwargs...)^c
# +(f::Function, c::Real) =
#     (args...; kwargs...) -> f(args...; kwargs...) + c
#     +(a, b) = +(b, a)`

δ(x,y) = ==(x,y)
export δ

struct ParameterList
    param_dict::Dict{Symbol, Float64}
end
export ParameterList

struct Operator{F}
    param::Symbol
    operator::F
    matrix::Matrix{ComplexF64}
end
export Operator

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
    for (i, bstate) ∈ enumerate(H.basis)
        for (j, bstate′) ∈ enumerate(H.basis)
            for p ∈ -1:1
                H.basis_tdms[i,j,p+2] = TDM(bstate, bstate′, p)
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
            end
        end
    end
    return nothing
end
export update_basis_tdms_m!

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
                H.tdms[j,i,p+2] += conj(H.tdms[i,j,p+2])
            end
        end
    end
    return nothing
end
export update_tdms!

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

function get_basis_tdms(basis::Vector{<:BasisState}, tdm_func::F) where F
    basis_tdms = zeros(ComplexF64, length(basis), length(basis), 3)
    for (i, bstate) ∈ enumerate(basis)
        for (j, bstate′) ∈ enumerate(basis)
            for p ∈ -1:1
                if i < j
                    basis_tdms[i,j,p+2] = tdm_func(bstate, bstate′, p)
                    basis_tdms[j,i,p+2] = conj(basis_tdms[i,j,p+2])
                elseif i == j
                    basis_tdms[i,j,p+2] = tdm_func(bstate, bstate′, p)
                end
            end
        end
    end
    return basis_tdms
end
export get_basis_tdms

function get_tdms_two_bases(basis::Vector{<:BasisState}, basis′::Vector{<:BasisState}, tdm_func::F) where F
    basis_tdms = zeros(ComplexF64, length(basis), length(basis′), 3)
    for (i, bstate) ∈ enumerate(basis)
        for (j, bstate′) ∈ enumerate(basis′)
            for p ∈ -1:1
                basis_tdms[i,j,p+2] = tdm_func(bstate, bstate′, p)
            end
        end
    end
    return basis_tdms
end
export get_tdms_two_bases

"""
    `states` and `states′` must share the same basis.
"""
function tdms_between_states!(d::Array{ComplexF64, 3}, basis_tdms::Array{ComplexF64, 3}, states::Vector{<:State}, states′::Vector{<:State})
    for i ∈ eachindex(states), j ∈ eachindex(states′)
        state  = states[i]
        state′ = states′[j]
        basis  = state.basis
        basis′ = state′.basis
        for p ∈ -1:1
            d[i,j,p+2] = 0.0
            for m ∈ eachindex(basis), n ∈ eachindex(basis′)
                d[i,j,p+2] += conj(state.coeffs[m]) * state′.coeffs[n] * basis_tdms[m,n,p+2]
            end
        end
    end
    return nothing
end
export tdms_between_states!

function add_to_H(H::Hamiltonian, param::Symbol, f::Function)
    operator = Expr(:call, :+, H.operator, f)
    matrix = matrix_from_operator(H.basis, f)
    param_dict = Dict(H.parameters.param_dict..., param => 0.0) # any new term added has its parameter value set to zero as default
    parameters = ParameterList(param_dict)
    operators = (; H.operators..., param => Operator(param, f, matrix))
    return Hamiltonian(H.basis, operator, parameters, H.states, operators, H.matrix, H.tdms, H.tdms_m, H.basis_tdms, H.basis_tdms_m)
end
export add_to_H

function +(H::Hamiltonian, O::Operator)
    H.operator = Expr(:call, :+, H.operator, O)
    H.parameters[O.param] = 0.0 # any new term added has its parameter value set to zero as default
    H.operators[O.param] = O.operator
    return H
end

function +(H::Hamiltonian, expr::Expr)
    param = expr.args[2]
    operator = eval(expr.args[3])
    matrix = matrix_from_operator(H.basis, operator)
    O = Operator(param, operator, matrix)
    return H + O
end
    
# function +(H1::Hamiltonian, H2::Hamiltonian)
#     Hamiltonian(H1.basis, )

macro params(arg)
    dict = Dict{Symbol, Float64}()
    for x in arg.args
        if x isa Expr
            param_name = x.args[1]
            param_value = x.args[2]
            dict[param_name] = eval(param_value)
        end
    end
    return ParameterList(dict)
end
export @params

import Base: getproperty
function getproperty(p::ParameterList, s::Symbol)
    param_dict = getfield(p, :param_dict)
    if s ∈ keys(param_dict)
        return getindex(param_dict, s)::Float64
    else
        return getfield(p, s)
    end
end
export getproperty

import Base: setproperty!
function setproperty!(p::ParameterList, s::Symbol, v::Float64)
    setindex!(getfield(p, :param_dict), v, s)
end
export setproperty!

macro make_operator(basis, operator_expr)
    return quote
        operators = NamedTuple()
        basis = $(esc(basis))
        # println($(length(operator_expr.args)))
        # println($(operator_expr.args[3]))
        for expr ∈ $(operator_expr.args)
            if expr isa Expr
                param = expr.args[2]
                println(param)
                operator = $(esc(operator_expr.args[2].args[3]))
                matrix = matrix_from_operator(basis, operator)
                Operator(param, operator, matrix)
                if param isa Symbol
                    println(operators)
                    operators = (param => 1,)
                    # operators = (; operators..., 1)
                end
            end
        end
        operators
    end
end
export make_operator

function unpack_operator(operator::Expr, basis::Vector{<:BasisState})
    operators = NamedTuple()
    if operator.args[1] != :+ # case if only one term is included in the Hamiltonian
        expr = operator
        param = expr.args[2]
        operator = eval(expr.args[3])
        matrix = matrix_from_operator(basis, operator)
        if param isa Symbol
            operators = (; operators..., param => Operator(param, operator, matrix))
        end
    else
        for expr ∈ operator.args
            if expr isa Expr
                param = expr.args[2]
                operator = eval(expr.args[3])
                matrix = matrix_from_operator(basis, operator)
                if param isa Symbol
                    operators = (; operators..., param => Operator(param, operator, matrix))
                end
            end
        end
    end
    return operators
end
export unpack_operator

function matrix_from_operator(basis, operator::F) where {F}
    m = zeros(ComplexF64, length(basis), length(basis))
    for i ∈ eachindex(basis), j ∈ eachindex(basis)
        m[i,j] = operator(basis[i], basis[j])
    end
    return m
end

# function evaluate_operator!(basis, o::Operator{F}) where F
#     for i ∈ eachindex(basis), j ∈ i:length(basis)
#         o.matrix[i,j] = o.operator(basis[i], basis[j])
#         o.matrix[j,i] = conj(o.matrix[i,j])
#     end
#     return nothing
# end
function evaluate_operator!(basis, o::Operator{F}) where F
    for i ∈ eachindex(basis), j ∈ eachindex(basis)
        o.matrix[i,j] = round( o.operator(basis[i], basis[j]), digits=20 )
    end
    return nothing
end
# export evaluate_operator!

@generated function evaluate_(H::Hamiltonian, _::Val{N}) where N
    quote
        Base.Cartesian.@nexprs $N i -> evaluate_operator!(H.basis, H.operators[i])
    end
end
export evaluate_

function add_operator_matrix!(H::Hamiltonian, operator::Operator)
    operator_val = getproperty(H.parameters, operator.param)
    for i ∈ eachindex(H.matrix)
        H.matrix[i] += operator_val * operator.matrix[i]
    end
    return nothing
end
export add_operator_matrix!

@generated function sum_operator_matrices!(H::Hamiltonian, _::Val{N}) where N
    quote
        Base.Cartesian.@nexprs $N i -> add_operator_matrix!(H, H.operators[i])
    end
end

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
export evaluate!
        
function solve!(H::Hamiltonian)
    es, vs = eigen(H.matrix)
    for i ∈ eachindex(H.states)
        H.states[i].E = real(es[i])
        H.states[i].coeffs = vs[:,i]
    end
    return nothing
end
export solve!

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

function findindex(QNs, QN::Symbol)
    exists = false
    i = 1
    for name in QNs
        if name == QN
            exists = true
            break
        else
            i += 1
        end
    end
    return (exists, i)
end

# Need to redefine so that this is type-stable...
function enumerate_states(state_type, QN_bounds1, QN_bounds2)
    basis_states1 = enumerate_states(state_type.types[1], QN_bounds1)
    basis_states2 = enumerate_states(state_type.types[2], QN_bounds2)
    basis = state_type[]
    for basis_state1 ∈ basis_states1
        for basis_state2 ∈ basis_states2
            push!(basis, state_type(basis_state1, basis_state2))
        end
    end
    return basis
end

function enumerate_states(state_type, QN_bounds)
    states = state_type[]
    QNs = [QN for QN ∈ fieldnames(state_type) if QN ∉ (:E, :constraints)]
    
    # Define a state with all QN = 0 to get the constraints for the QNs
    η = NamedTuple([QN => 0 for QN in QNs])
    
    enumerate_states(η, states, state_type, QNs, QN_bounds, 1)

    return states
end

function enumerate_states(η, states, state_type, QNs, QN_bounds, idx, max_states=1000)

    if length(states) > max_states
        return false
    end

    if idx == length(QNs) + 1
        new_state = state_type(; η...)
        push!(states, new_state)
        return nothing
    end
    
    # Check if quantum number has been given bounds; else apply constraints
    iterated_QN = QNs[idx]
    QN_constraints = state_type(; η...).constraints

    if iterated_QN ∈ keys(QN_constraints)
        QN_constraint = QN_constraints[iterated_QN]
        QN_constraint_bounds = eval(QN_constraint)
    end
    
    if iterated_QN ∈ keys(QN_bounds)
        bounds = QN_bounds[iterated_QN]
        for i ∈ bounds
            if iterated_QN ∉ keys(QN_constraints) || (i ∈ QN_constraint_bounds)
                η′ = (; η..., iterated_QN => i)
                enumerate_states(η′, states, state_type, QNs, QN_bounds, idx + 1)
            end
        end
    elseif iterated_QN ∈ keys(QN_constraints)
        for i ∈ QN_constraint_bounds
            η′ = (; η..., iterated_QN => i)
            enumerate_states(η′, states, state_type, QNs, QN_bounds, idx + 1)
        end
    else
        enumerate_states(η, states, state_type, QNs, QN_bounds, idx + 1)
    end
    
    return nothing
end
export enumerate_states

# function enumerate_states(η, QNs, states, state_type, QN_bounds, idx, max_states=1000)

#     if length(states) > max_states
#         return false
#     end

#     return_val = true
#     try
#         new_state = state_type(; η...)
#         return_val = true
#     catch
#         return_val = false
#     end

#     if idx == length(QNs) + 1
#         # new_state = state_type(; η...)
#         try
#             new_state = state_type(; η...)
#             push!(states, new_state)
#             return true
#         catch
#             return false
#         end
#     end
    
#     # Check if quantum number has been given bounds
#     iterated_QN = QNs[idx]
#     QN_bounded = false
#     if iterated_QN in keys(QN_bounds)
#         QN_bounded = true
#         bounds = QN_bounds[iterated_QN]
#         if η[iterated_QN] ∉ bounds
#             η = (; η..., iterated_QN => bounds[1])
#         end
#     end
            
#     i = η[iterated_QN] + 1
#     keep_iterating = true
#     # while keep_iterating && ((!QN_bounded) || (QN_bounds[iterated_QN][1] <= i <= QN_bounds[iterated_QN][2]))
#     while keep_iterating && ((!QN_bounded) || (i ∈ bounds))
#         η′ = (; η..., iterated_QN => i)
# #         println(η′)
# #         println(iterated_QN, i)
#         keep_iterating = enumerate_states(η′, QNs, states, state_type, QN_bounds, idx + 1, max_states)
#         i += 1
#     end
    
#     i = η[iterated_QN] - 1
#     keep_iterating = true
#     # while keep_iterating && ((!QN_bounded) || (QN_bounds[iterated_QN][1] <= i <= QN_bounds[iterated_QN][2]))
#     while keep_iterating && ((!QN_bounded) || (i ∈ bounds))
#         η′ = (; η..., iterated_QN => i)
# #         println(η′)
# #         println(iterated_QN, i)
#         keep_iterating = enumerate_states(η′, QNs, states, state_type, QN_bounds, idx + 1, max_states)
#         i -= 1
#     end
    
#     return return_val
# end
# export enumerate_states

import Base.==
function ==(state::State, state′::State)
    state_type = typeof(state)
    for field in fieldnames(state_type)
        if getfield(state, field) != getfield(state′, field)
            return false
        end
    end
    return true
end

function ==(state::BasisState, state′::BasisState)
    state_type = typeof(state)
    for field in fieldnames(state_type)
        if getfield(state, field) != getfield(state′, field)
            return false
        end
    end
    return true
end

export ==

function DiagonalOperator(state::BasisState, state′::BasisState)
    if state == state′
        return 1.0
    else
        return 0.0
    end
end
export DiagonalOperator

mutable struct Transition
    ground_state::State
    excited_state::State
    frequency::Float64
    tdm::ComplexF64
end
export Transition

function compute_transitions(H::Hamiltonian, p, threshold=1e-8)
    transitions = Transition[]
    
    for (i, basis_state) ∈ enumerate(H.basis)
        for (j, basis_state′) ∈ enumerate(H.basis)
            H.tdms[i,j,p] = TDM(basis_state, basis_state′, p)
        end
    end
            
    for (i, state) ∈ enumerate(H.states)
        for (j, state′) ∈ enumerate(H.states)
            if state′.E > state.E
                tdm = state.coeffs ⋅ (H.tdms * state′.coeffs)
                if norm(tdm) > threshold
                    transition = Transition(state, state′, state′.E - state.E, tdm)
                    push!(transitions, transition)
                    H.tdms
                end
            end
        end
    end
    return transitions
end

"""
    Compute the transitions between two sets of states, `states` and `states′`.
"""
function compute_transitions(states::Vector{<:State}, states′::Vector{<:State}, p; threshold=1e-8, compute_tdms=true)
    transitions = Transition[]
    basis = states[1].basis

    tdms = zeros(length(basis), length(basis))
    if compute_tdms
        for (i, basis_state) ∈ enumerate(basis)
            for (j, basis_state′) ∈ enumerate(basis)
                tdms[i,j] = TDM(basis_state, basis_state′, p)
            end
        end
    end

    for state ∈ states
        for state′ ∈ states′
            if state′.E > state.E
                tdm = state.coeffs ⋅ (tdms * state′.coeffs)
                f = state′.E - state.E
                if (norm(tdm) > threshold || ~compute_tdms) && abs(f) > 1
                    transition = Transition(state, state′, f, tdm)
                    push!(transitions, transition)
                end
            end
        end
    end
    return transitions
end
export compute_transitions

ground_state(transition::Transition) = transition.ground_state
excited_state(transition::Transition) = transition.excited_state
frequency(transition::Transition) = transition.frequency
tdm(transition::Transition) = transition.tdm
export ground_state
export excited_state

function transitions_table(transitions::Vector{Transition}, relabelling_states=nothing::Union{Nothing, Vector{<:State}}; threshold=1e-8)
    
    ground_states = ground_state.(transitions)
    excited_states = excited_state.(transitions)
    frequencies = frequency.(transitions)
    tdms = tdm.(transitions)
    
    if ~isnothing(relabelling_states)
        relabelled_ground_states = relabelling_states[getfield.(ground_states, :idx)]
        relabelled_excited_states = relabelling_states[getfield.(excited_states, :idx)]
    else
        relabelled_ground_states = ground_states
        relabelled_excited_states = excited_states
    end
    
    df1 = states_table(relabelled_ground_states, threshold=1e-1, dominant_state_only=true)
    df2 = states_table(relabelled_excited_states, threshold=1e-1, dominant_state_only=true)
    
    df_transitions = hcat(df1, df2, makeunique=true)
    df_transitions[!, :f] = frequencies
    df_transitions[!, :tdm] = real.(tdms)
    sort!(df_transitions, :f)
    
    return df_transitions
end
export transitions_table

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

function calculate_state_overlaps!(states, states′, overlaps)
    for i ∈ eachindex(states)
        for j ∈ eachindex(states′)
            state = states[i]
            state′ = states′[j]
            overlaps[i,j] = norm(state ⋅ state′)
        end
    end
    return nothing
end  

function calculate_state_overlaps(states, states′)
    overlaps = zeros(Float64, length(states), length(states))
    calculate_state_overlaps!(states, states′, overlaps)
    return overlaps
end
export calculate_state_overlaps

function tracked_idxs!(overlaps, states, prev_states, tracked_idxs)
    for i ∈ axes(overlaps, 2)
        overlaps_row = @view overlaps[:,i] # overlap for previous state `i` with all states of current iteration 
        max_idx = argmax(overlaps_row)
        tracked_idxs[i] = max_idx
        overlaps[max_idx,:] .= 0.0
        # found_idx = false
        # j = 1
        # while !found_idx && j <= size(overlaps,1)
        #     max_idx = argmax(overlaps_row)
        #     # println(abs(energy(prev_states[i]) - energy(states[tracked_idxs[i]])), " ", abs(energy(prev_states[i]) - energy(states[max_idx])))
        #     # println(energy(prev_states[i]) - energy(states[tracked_idxs[i]]))
        #     # println(energy(prev_states[i]), " ", energy(prev_states[i]) - energy(states[max_idx]), " ", j)
        #     # println(" ")
        #     # Don't allow the tracked index to change if the change is energy is too large
        #     # (This is meant to avoid discrete jumps between states far from each other that may look similar close to an avoided crossing.)
        #     if abs(energy(prev_states[i]) - energy(states[max_idx])) < (10 * abs(energy(prev_states[i]) - energy(states[i])))
        #         tracked_idxs[i] = max_idx
        #         overlaps[max_idx,:] .= 0.0
        #         found_idx = true
        #     else
        #         overlaps[max_idx,i] = 0.0
        #     end
        #     j += 1
        # end
        # println(abs(energy(prev_states[i]) - energy(states[i])), " ", abs(energy(prev_states[i]) - energy(states[tracked_idxs[i]])))
        # if !found_idx
            # println(found_idx, i)
        # end
    end
    return nothing
end
export tracked_idxs!

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
# export scan_single_parameter
# export scan_parameters

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