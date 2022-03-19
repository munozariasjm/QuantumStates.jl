module QuantumStates

include("States.jl")
include("HundsCaseA.jl")
include("HundsCaseB.jl")
include("Printing.jl")
include("Operators.jl")

using HalfIntegers
using Parameters
using NamedTupleTools
using LinearAlgebra

function states_from_basis(basis::Vector{<:BasisState})::Vector{<:State}
    """
    Returns a vector of states with a one-to-one correspondence with the basis.
    """
    states = State[]
    for i in eachindex(basis)
        coeffs = zeros(ComplexF64, length(basis))
        coeffs[i] = 1
        push!(states, State(basis, coeffs))
    end
    return states
end
export states_from_basis

δ(x,y) = ==(x,y)
export δ

# @with_kw mutable struct States
#     basis::Array{State, 1}
#     states::Array{ComplexF64, 2} = diagm(ones(Float64, length(basis)))
# end
# export States

@with_kw mutable struct Hamiltonian{T}
    basis::Vector{<:BasisState}
    H_operator::T
    M::Array{ComplexF64, 2} = [H_operator(state, state′) for state in basis, state′ in basis]
#     M::Array{ComplexF64, 2} = zeros(ComplexF64, length(basis), length(basis))
end
export Hamiltonian

function solve(H::Hamiltonian; kwargs...)
    basis = H.basis
    H_basis = zeros(ComplexF64, length(basis), length(basis))

    args = [kwarg.second for kwarg in kwargs]
    if !isempty(args)
        H_operator = H.H_operator(args...)
    else
        H_operator = H.H_operator
    end
    for (i, state) in enumerate(basis)
        for (j, state′) in enumerate(basis)
            H_basis[i,j] = H_operator(state, state′)
        end
    end
    H.M = H_basis #states * H_basis * states'
#     eigvals, eigvecs = eigen(H.M)
    return eigen(H.M)
end
export solve

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

function enumerate_states(η, state_type, QN_bounds)
    states = state_type[]
    QNs = fieldnames(state_type)

    bounds_exist = true
    earliest_bound_idx = length(QNs)
    for bound_QN in keys(QN_bounds)
        bound_QN_exists, QN_idx_in_state = findindex(QNs, bound_QN)
        if bound_QN_exists
            if QN_idx_in_state < earliest_bound_idx
                earliest_bound_idx = QN_idx_in_state
            end
        else
            bounds_exist = false
        end
    end

    if bounds_exist
        enumerate_states(η, QNs, states, state_type, QN_bounds, earliest_bound_idx)
    else
        error("Bounded quantum number does not exist in this basis.")
    end
    return states
end

function enumerate_states(η, QNs, states, state_type, QN_bounds, idx, max_states=1000)

    if length(states) > max_states
        return false
    end

    return_val = true
    try
        new_state = state_type(; η...)
        return_val = true
    catch
        return_val = false
    end

    if idx == length(QNs) + 1
        # new_state = state_type(; η...)
        try
            new_state = state_type(; η...)
            push!(states, new_state)
            return true
        catch
            return false
        end
    end

    # Check if quantum number has been given bounds
    iterated_QN = QNs[idx]
    QN_bounded = false
    if iterated_QN in keys(QN_bounds)
        QN_bounded = true
        bounds = QN_bounds[iterated_QN]
        if η[iterated_QN] ∉ bounds
            η = (; η..., iterated_QN => bounds[1])
        end
    end

    keep_iterating = enumerate_states(η, QNs, states, state_type, QN_bounds, idx + 1, max_states)

    i = η[iterated_QN] + 1
    keep_iterating = true
    # while keep_iterating && ((!QN_bounded) || (QN_bounds[iterated_QN][1] <= i <= QN_bounds[iterated_QN][2]))
    while keep_iterating && ((!QN_bounded) || (i ∈ bounds))
        η′ = (; η..., iterated_QN => i)
#         println(η′)
#         println(iterated_QN, i)
        keep_iterating = enumerate_states(η′, QNs, states, state_type, QN_bounds, idx + 1, max_states)
        i += 1
    end

    i = η[iterated_QN] - 1
    keep_iterating = true
    # while keep_iterating && ((!QN_bounded) || (QN_bounds[iterated_QN][1] <= i <= QN_bounds[iterated_QN][2]))
    while keep_iterating && ((!QN_bounded) || (i ∈ bounds))
        η′ = (; η..., iterated_QN => i)
#         println(η′)
#         println(iterated_QN, i)
        keep_iterating = enumerate_states(η′, QNs, states, state_type, QN_bounds, idx + 1, max_states)
        i -= 1
    end

    return return_val
end
export enumerate_states

# function convert_basis(H::Hamiltonian, basis′, overlap)

#     P = deepcopy(H.M)
#     P .= 0.0

#     basis = H.basis
#     for (i, state) in enumerate(basis)
#         for (j, state′) in enumerate(basis′)
#             P[i,j] = overlap(state, state′)
#         end
#     end
#     H.M .= P' * H.M * P
#     H.basis = basis′

#     return P
# end
# function convert_basis(states::States, basis′, overlap)

#     new_states = deepcopy(states)
#     new_states.basis = basis′

#     basis = states.basis
#     states = states.states

#     # Calculate change-of-basis matrix
#     P = zeros(Float64, length(basis′), length(basis))
#     for (i, state) in enumerate(basis′)
#         for (j, state′) in enumerate(basis)
#             P[i,j] = overlap(state, state′)
#         end
#     end

#     new_states.states = round.(P * states, digits=10)

#     return new_states
# end
# export convert_basis

function get_subspace(states, QN_bounds)
    new_basis = []

    QN_bound_symbol = keys(QN_bounds)[1]
    bounds = QN_bounds[1][1], QN_bounds[1][2]

    for (i, state) in enumerate(states.basis)
        state_val = getfield(state, QN_bound_symbol)
        if bounds[1] <= state_val <= bounds[2]
            push!(states_bare, state)
        end
    end

    return States(basis=states_bare)
end
export get_subspace

@with_kw mutable struct HundsCaseA_Parity <: BasisState
    S::HalfInteger
    I::HalfInteger
    Λ::HalfInteger
    Σ::HalfInteger
    Ω::HalfInteger
    P::HalfInteger
    J::HalfInteger
    F::HalfInteger
    M::HalfInteger
    function HundsCaseA_Parity(S, I, Λ, Σ, Ω, P, J, F, M)
        if abs(Σ) > S
            error("|Σ| > S")
        elseif !(abs(I - J) <= F <= abs(J + I))
            error("F < |I - J| or F > |I + J|")
#         elseif !(abs(Ω) == abs(Σ + Λ))
#             error("|Ω| != |Σ + Λ|")
        elseif abs(Ω) > J
            error("|Ω| > J")
        elseif abs(M) > F
            error("|M| > F")
        elseif abs(P) != 1
            error("P != ±1")
        elseif Λ < 0 || Σ < 0 || Ω < 0
            error("Λ, Σ, Ω must be absolute values.")
        end
        return new(S, I, Λ, Σ, Ω, P, J, F, M)
    end
end
export HundsCaseA_Parity

# Extend the base "+" function such that we can simply add matrix elements before calling them
import Base.+, Base.*

+(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) + g(args...; kwargs...)
*(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) * g(args...; kwargs...)
*(c::Float64, f::Function) =
    (args...; kwargs...) -> c * f(args...; kwargs...)
export +, *

function overlap(state::HundsCaseB, state′::HundsCaseA)
    # Eq. (6.149) in Brown & Carrington
    S, I, Λ, N, J, F, M = state.S, state.I, state.Λ, state.N, state.J, state.F, state.M
    S′, I′, Λ′, Σ′, Ω′, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    return (-1)^(N - S + Ω′) * sqrt(2N + 1) * wigner3j_(J, S, N, Ω′, -Σ′, -Λ) * δ(J, J′) * δ(M, M′) * δ(F, F′)
end
overlap(state::HundsCaseA, state′::HundsCaseB) = overlap(state′, state)

function overlap(state::HundsCaseA, state′::HundsCaseA_Parity)
    S, I, Λ, Σ, Ω, J, F, M =
        state.S, state.I, state.Λ, state.Σ, state.Ω, state.J, state.F, state.M
    S′, I′, Λ′, Σ′, Ω′, P′, J′, F′, M′ =
        state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.P, state′.J, state′.F, state′.M

    if Λ < 0
        parity_factor = P′ * (-1)^(J - S)
    else
        parity_factor = 1
    end

    return (1/√2) * parity_factor * δ(S, S′) * δ(I, I′) * δ(J, J′) * δ(F, F′) * δ(M, M′) *
        δ(abs(Ω), Ω′) * δ(abs(Σ), Σ′) * δ(abs(Λ), Λ′)
end
overlap(state::HundsCaseA_Parity, state′::HundsCaseA) = overlap(state′, state)
export overlap

end
