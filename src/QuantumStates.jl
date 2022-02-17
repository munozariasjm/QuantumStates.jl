module QuantumStates

using HalfIntegers 
using Parameters
using NamedTupleTools
using LinearAlgebra
import WignerSymbols: wigner3j, wigner6j

δ(x,y) = ==(x,y)
export δ

abstract type State end
export State

function print_nice(state::State)
    for field in fieldnames(typeof(state))
        println(field, ": ", getfield(state, field))
    end
    return nothing
end
export print_nice

@with_kw mutable struct Hamiltonian
    basis::Array{State, 1}
    M::Array{ComplexF64, 2} = zeros(ComplexF64, length(basis), length(basis))
end
export Hamiltonian

@with_kw mutable struct States
    basis::Array{State, 1}
    states::Array{Float64, 2} = diagm(ones(Float64, length(basis)))
end
export States

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

    QN_bound_symbol = keys(QN_bounds)[1]
    bounds = QN_bounds[1][1], QN_bounds[1][2]

    (bound_QN_exists, bound_QN_idx) = findindex(QNs, QN_bound_symbol)

    η = @eval ($η..., $QN_bound_symbol = $bounds[1])

    if bound_QN_exists
        enumerate_states(η, QNs, states, state_type, bound_QN_idx, bounds, bound_QN_idx)
    else
        error("Bounded quantum number does not exist in this basis.")
    end
    return states
end

function enumerate_states(η, QNs, states, state_type, bound_QN_idx, bounds, idx=bound_QN_idx, max_states=1000)

    if length(states) > max_states
        return false
    end

    if idx == length(fieldnames(state_type)) + 1
        # new_state = state_type(; η...)
        try
            new_state = state_type(; η...)
            push!(states, new_state)
            return true
        catch
            return true
        end
    end

    prev_QN = QNs[idx-1]
    iterated_QN = QNs[idx]

    if idx == bound_QN_idx
        lower_limit = bounds[1]
        upper_limit = bounds[2]
    else
        upper_limit = abs(η[prev_QN]) + 1
        lower_limit = -upper_limit
    end

    i = η[iterated_QN]
    keep_iterating = true
    while keep_iterating && i <= upper_limit
        #println(η)
        η′ = (; η..., iterated_QN => i)
        keep_iterating = enumerate_states(η′, QNs, states, state_type, bound_QN_idx, bounds, idx + 1, max_states)
        i += 1
    end

    i = η[iterated_QN] - 1
    keep_iterating = true
    while keep_iterating && lower_limit <= i
        η′ = (; η..., iterated_QN => i)
        keep_iterating = enumerate_states(η′, QNs, states, state_type, bound_QN_idx, bounds, idx + 1, max_states)
        i -= 1
    end

#     i = η[iterated_QN]
#     keep_iterating = true
#     while keep_iterating #|| (idx == bound_QN_idx && bounds[1] <= η[QNs[bound_QN_idx]] <= bounds[2])
#         η′ = (; η..., iterated_QN => i)
#         keep_iterating = enumerate_states(η′, QNs, states, state_type, bound_QN_idx, bounds, idx + 1, max_states)
#         i += 1
#     end

#     i = η[iterated_QN] - 1
#     keep_iterating = true
#     while keep_iterating #|| (idx == bound_QN_idx && bounds[1] <= η[QNs[bound_QN_idx]] <= bounds[2])
#         η′ = (; η..., iterated_QN => i)
#         keep_iterating = enumerate_states(η′, QNs, states, state_type, bound_QN_idx, bounds, idx + 1, max_states)
#         i -= 1
#     end

    return true
end
export enumerate_states

function convert_basis(H::Hamiltonian, basis′, overlap)

    P = deepcopy(H.M)
    P .= 0.0

    basis = H.basis
    for (i, state) in enumerate(basis)
        for (j, state′) in enumerate(basis′)
            P[i,j] = overlap(state, state′)
        end
    end
    H.M .= P' * H.M * P
    H.basis = basis′

    return P
end
function convert_basis(states::States, basis′, overlap)

    new_states = deepcopy(states)
    new_states.basis = basis′
    
    basis = states.basis
    states = states.states
    
    # Calculate change-of-basis matrix
    P = zeros(Float64, length(basis′), length(basis))
    for (i, state) in enumerate(basis′)
        for (j, state′) in enumerate(basis)
            P[i,j] = overlap(state, state′)
        end
    end
    
    new_states.states = round.(P * states, digits=10)
    
    return new_states
end
export convert_basis

function get_subspace(H, QN_bounds)

    idxs = []
    
    QN_bound_symbol = keys(QN_bounds)[1]
    bounds = QN_bounds[1][1], QN_bounds[1][2]
    
    for (i, state) in enumerate(H.basis)
        state_val = getfield(state, QN_bound_symbol)
        if bounds[1] <= state_val <= bounds[2]
            push!(idxs, i)
        end
    end
    
    H_subspace = deepcopy(H)
    H_subspace.M = H.M[idxs, idxs]
    H_subspace.basis = H.basis[idxs]
    
    return H_subspace
end
export get_subspace

@with_kw mutable struct HundsCaseB <: State
    S::HalfInteger 
    I::HalfInteger
    Λ::HalfInteger 
    N::HalfInteger 
    J::HalfInteger 
    F::HalfInteger
    M::HalfInteger
    function HundsCaseB(S, I, Λ, N, J, F, M)
        if abs(Λ) > N
            error("|Λ| > N")
        elseif !(abs(N - 1//2) <= J <= N + 1//2)
            error("J > N + S")
        elseif !(abs(J - 1//2) <= F <= J + 1//2)
            error("F > I + J")
        elseif abs(M) > F
            error("|M| > F")
        end
        return new(S, I, Λ, N, J, F, M)
    end
end
export HundsCaseB

@with_kw mutable struct HundsCaseA <: State
    S::HalfInteger
    I::HalfInteger
    Λ::HalfInteger
    Σ::HalfInteger
    Ω::HalfInteger
    J::HalfInteger
    F::HalfInteger
    M::HalfInteger
    function HundsCaseA(S, I, Λ, Σ, Ω, J, F, M)
        if abs(Σ) > S
            error("|Σ| > S")
        elseif !(abs(I - J) <= F <= abs(J + I))
            error("F < |I - J| or F > |I + J|")
        elseif !(Ω == Σ + Λ)
            error("Ω != Σ + Λ")
        elseif abs(Ω) > J
            error("|Ω| > J")
        elseif abs(M) > F
            error("|M| > F")
        end
        return new(S, I, Λ, Σ, Ω, J, F, M)
    end
end
export HundsCaseA

@with_kw mutable struct HundsCaseA_Parity <: State
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

# Define matrix elements for Hund's case (a)
function Rotation(state::HundsCaseA, state′::HundsCaseA)
    S, I, Λ, Σ, Ω, J, F, M = state.S, state.I, state.Λ, state.Σ, state.Ω, state.J, state.F, state.M
    S′, I′, Λ′, Σ′, Ω′, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    return (J * (J + 1) + S * (S + 1) - Ω^2 - Σ^2) * 
        δ(Λ, Λ′) * δ(Ω, Ω′) * δ(Σ, Σ′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
end

function SpinOrbit(state::HundsCaseA, state′::HundsCaseA)
    S, I, Λ, Σ, Ω, J, F, M = state.S, state.I, state.Λ, state.Σ, state.Ω, state.J, state.F, state.M
    S′, I′, Λ′, Σ′, Ω′, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    return Λ * Σ * δ(Λ, Λ′) * δ(Ω, Ω′) * δ(Σ, Σ′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
end
export SpinOrbit

function Hyperfine_IL(state::HundsCaseA, state′::HundsCaseA)
    S, I, Λ, Σ, Ω, J, F, M = state.S, state.I, state.Λ, state.Σ, state.Ω, state.J, state.F, state.M
    S′, I′, Λ′, Σ′, Ω′, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    return Λ * (-1)^(J′ + I + F + J - Ω) * 
        sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1)) *
        wigner6j(J′, I, F, I, J, 1) * 
        wigner3j(J, 1, J′, -Ω, 0, Ω′) * 
        δ(Λ, Λ′) * δ(Σ, Σ′) * δ(F, F′) * δ(M, M′)
end
export Hyperfine_IL

function Hyperfine_IF(state::HundsCaseA, state′::HundsCaseA)
    S, I, Λ, Σ, Ω, J, F, M = state.S, state.I, state.Λ, state.Σ, state.Ω, state.J, state.F, state.M
    S′, I′, Λ′, Σ′, Ω′, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    return sum(-1)^(I + J′ + F + S - Σ + J - Ω) * 
        sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1)) *
        wigner6j(J′, I, F, I, J, 1) *
        δ(Λ, Λ′) * δ(F, F′) * δ(M, M′) *
        sum(
            wigner3j(J, 1, J′, -Ω, q, Ω′) *
            wigner3j(S, 1, S, -Σ, q, Σ′)
            for q in -1:1
        )
end
export Hyperfine_IF

# Define matrix elements for Hund's case (b)
function unpack(state::HundsCaseB)
    return (state.S, state.I, state.N, state.Λ, state.J, state.F, state.M)
end
export unpack

function Rotation(state::HundsCaseB, state′::HundsCaseB)
    S, I, N, Λ, J, F, M = unpack(state)
    S′, I′, N′, Λ′, J′, F′, M′ = unpack(state′)
    return N * (N + 1) * δ(Λ, Λ′) * δ(N, N′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
end
export Rotation

function SpinRotation(state::HundsCaseB, state′::HundsCaseB)
    S, I, N, Λ, J, F, M = unpack(state)   
    S′, I′, N′, Λ′, J′, F′, M′ = unpack(state′)
    return (-1)^(N + J + S) * sqrt(S * (S + 1) * (2S + 1)) * sqrt(N * (N + 1) * (2N + 1)) * 
        wigner6j(S, N, J, N, S, 1) *
        δ(Λ, Λ′) * δ(N, N′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
end
export SpinRotation

function Hyperfine_IS(state::HundsCaseB, state′::HundsCaseB)
    S, I, N, Λ, J, F, M = unpack(state)   
    S′, I′, N′, Λ′, J′, F′, M′ = unpack(state′)
    return (-1)^(J′ + F + I + J + N + S + 1) *
        sqrt( (2J′ + 1) * (2J + 1) * S * (S + 1) * (2S + 1) * I * (I + 1) * (2I + 1) ) *
        wigner6j(I, J′, F, J, I, 1) *
        wigner6j(J, S, N, S, J′, 1) *
        δ(Λ, Λ′) * δ(N, N′) * δ(F, F′) * δ(M, M′)
end
export Hyperfine_IS

function Hyperfine_IK(state::HundsCaseB, state′::HundsCaseB)
    S, I, N, Λ, J, F, M = unpack(state)   
    S′, I′, N′, Λ′, J′, F′, M′ = unpack(state′)
    return (-1)^(F + I + N′ + S + 2J + 1 + N - Λ) * sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * (2N + 1) * (2N′ + 1)) *
        wigner6j(I, J, F, J′, I, 1) * 
        wigner6j(N, J, S, J′, N′, 1) * 
        wigner3j(N′, 1, N, -Λ, 0, Λ) *
        δ(Λ, Λ′) * δ(F, F′) * δ(M, M′)
end
export Hyperfine_IK

function Hyperfine_SK(state::HundsCaseB, state′::HundsCaseB)
    S, I, N, Λ, J, F, M = unpack(state)   
    S′, I′, N′, Λ′, J′, F′, M′ = unpack(state′)
    return (-1)^(2N + J + S - Λ) * sqrt(S * (S + 1) * (2S + 1) * (2N + 1) * (2N′ + 1)) *
        wigner6j(S, N, J, N′, S, 1) * 
        wigner3j(N′, 1, N, -Λ, 0, Λ) *
        δ(Λ, Λ′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
end
export Hyperfine_SK

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
    return (-1)^(N - S + Ω′) * sqrt(2N + 1) * wigner3j(J, S, N, Ω′, -Σ′, -Λ) * δ(J, J′) * δ(M, M′) * δ(F, F′)
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
