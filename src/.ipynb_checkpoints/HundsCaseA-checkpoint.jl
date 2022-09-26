using Parameters

abstract type HundsCaseA <: BasisState end
export HundsCaseA

@composite Base.@kwdef struct HundsCaseA_Rot <: HundsCaseA
    E::Float64 = 0.0
    I::HalfInt
    S::HalfInt
    Σ::HalfInt
    Λ::HalfInt
    J::HalfInt
    Ω::HalfInt
    F::HalfInt
    M::HalfInt
    constraints = (
        Σ = -S:S,
        Ω = max(-J, Λ + Σ):min(J, Λ + Σ),
        F = abs(J - I):abs(J + I),
        M = -F:F
    )
end
export HundsCaseA_Rot

function unpack(state::HundsCaseA)
    return (state.I, state.S, state.Λ, state.J, state.Ω, state.Σ, state.F, state.M)
end
export unpack

function Rotation(state::HundsCaseA, state′::HundsCaseA)
    I,  S,  Λ,  J,  Ω,  Σ,  F,  M  = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return (J * (J + 1) + S * (S + 1) + Λ^2 - 2 * Ω^2 + 2Λ * Σ) *
        δ(Λ, Λ′) * δ(Σ, Σ′) * δ(Λ, Λ′) * δ(J, J′) * δ(Ω, Ω′) * δ(F, F′) * δ(M, M′)
end

function SpinOrbit(state::HundsCaseA, state′::HundsCaseA)
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return Λ * Σ * 
        δ(Λ, Λ′) * δ(Σ, Σ′) * δ(Λ, Λ′) * δ(J, J′) * δ(Ω, Ω′) * δ(F, F′) * δ(M, M′)
end
export SpinOrbit

function SpinUncoupling(state::HundsCaseA, state′::HundsCaseA)
    # Brown & Carrington, eq. (8.364)
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return -2 * (-1)^(J - Ω + S - Σ) * 
        sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
        sum(
            wigner3j_(J, 1, J, -Ω, q, Ω′) * 
            wigner3j_(S, 1, S, -Σ, q, Σ′)
            for q in [-1,1]
        ) *
        δ(Λ, Λ′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
end
export SpinUncoupling

function ΛDoubling_q(state::HundsCaseA, state′::HundsCaseA)
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return (-1)^(J - Ω) *
        (1 / (2 * sqrt(6))) *
        sqrt( (2J - 1) * (2J) * (2J + 1) * (2J + 2) * (2J + 3) ) *
        sum(
            wigner3j_(J, 2, J, -Ω, -2q, Ω′)
            for q in [-1,1]
        ) *
        δ(abs(Λ′ - Λ), 2) *
        δ(J, J′) * δ(F, F′) * δ(M, M′)
end
export ΛDoubling_q

function ΛDoubling_p2q(state::HundsCaseA, state′::HundsCaseA)
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return (-1)^(J - Ω + S - Σ) *
        sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
        sum(
            wigner3j_(J, 1, J, -Ω, -q, Ω′) *
            wigner3j_(S, 1, S, -Σ, q, Σ′)
            for q in [-1,1]
        ) *
        δ(abs(Λ′ - Λ), 2) *
        δ(J, J′) * δ(F, F′) * δ(M, M′)
    end
export ΛDoubling_p2q
        
function Hyperfine_IL(state::HundsCaseA, state′::HundsCaseA)
    # Brown and Carrington, eq. (8.372)
    # Orbital hyperfine interaction
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return Λ * (-1)^(J′ + I + F + J - Ω) * 
        sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1)) *
        wigner6j_(I, J′, F, J, I, 1) * 
        wigner3j_(J, 1, J′, -Ω, 0, Ω′) * 
        δ(Λ, Λ′) * δ(F, F′) * δ(M, M′)
end
export Hyperfine_IL

function Hyperfine_IF(state::HundsCaseA, state′::HundsCaseA)
    # Fermi contact interaction
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return (-1)^(I + J′ + F + S - Σ + J - Ω) * 
        sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1)) *
        wigner6j_(I, J′, F, J, I, 1) *
        sum(
            wigner3j_(J, 1, J′, -Ω, q, Ω′) *
            wigner3j_(S, 1, S, -Σ, q, Σ′)
            for q in -1:1
        ) *
        δ(F, F′) * δ(M, M′)
end
export Hyperfine_IF
    
function Hyperfine_Dipolar_c(state::HundsCaseA, state′::HundsCaseA)
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return sqrt(30) * (1/3) * (-1)^(I + J′ + F) * (-1)^(J - Ω) * (-1)^(S - Σ) *
        wigner6j_(I, J′, F, J, I, 1) *
        sqrt( I * (I + 1) * (2I + 1) ) *
        sqrt( (2J + 1) * (2J′ + 1) ) *
        sqrt( S * (S + 1) * (2S + 1) ) *
        sum(
            (-1)^q * 
            wigner3j_(J, 1, J′, -Ω, q, Ω) *
            sum(
                wigner3j_(1, 2, 1, q′, 0, -q) *
                wigner3j_(S, 1, S, -Σ, q′, Σ′)
                for q′ in -1:1
            ) for q in -1:1
        ) * 
        δ(F, F′) * δ(M, M′)
end
export Hyperfine_Dipolar_c

function Hyperfine_Dipolar_d(state::HundsCaseA, state′::HundsCaseA)
    I, S, Λ, J, Ω, Σ, F, M = unpack(state)
    I′, I′, Λ′, J′, Ω′, Σ′, F′, M′ = unpack(state′)
    return sqrt(30) * (1/2) * sqrt(3/2) * (2/3) * (-1)^(I + J′ + F) * (-1)^(J - Ω) * (-1)^(S - Σ) *
        wigner6j_(I, J′, F, J, I, 1) *
        sqrt( I * (I + 1) * (2I + 1) ) *
        sqrt( (2J + 1) * (2J′ + 1) ) *
        sqrt( S * (S + 1) * (2S + 1) ) *
        sum(
            (-1)^q * 
            wigner3j_(J, 1, J′, -Ω, q, Ω) *
            sum(
                (
                    wigner3j_(1, 2, 1, q′, 2, -q) +
                    wigner3j_(1, 2, 1, q′, -2, -q)
                ) *
                wigner3j_(S, 1, S, -Σ, q′, Σ′)
                for q′ in -1:1
            ) for q in -1:1
        ) * 
        δ(F, F′) * δ(M, M′)
end
export Hyperfine_Dipolar_d

function TDM(state::HundsCaseA, state′::HundsCaseA, p::Integer)
    S, I, Λ, Σ, Ω, J, F, M = state.S, state.I, state.Λ, state.Σ, state.Ω, state.J, state.F, state.M
    S′, I′, Λ′, Σ′, Ω′, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    
    return (-1)^(F′ - M′) * wigner3j(F′, 1, F, -M′, p, M) * (-1)^(F + J′ + I + 1) * sqrt((2F′ + 1) * (2F + 1)) *
        δ(Σ, Σ′) * (-1)^(J′ - Ω′) * wigner6j(J, F, I, F′, J′, 1) * sqrt((2J′ + 1) * (2J + 1)) *
        sum(
            wigner3j(J′, 1, J, -Ω′, q, Ω) for q in -1:1
        )
end
TDM(state::HundsCaseA, state′::HundsCaseA) = sum(TDM(state, state′, p) for p ∈ -1:1)
export TDM

# function conversion(state::HundsCaseA, basis_type)
#     """
#     Returns a tuple of state vector of basis states
#     """
#     new_basis_states = HundsCaseB[]
#     coeffs = ComplexF64[]
    
#     for Σ in -S:S
#         Ω = Σ + Λ
#         coeff = (-1)^(J - S + Λ) * sqrt(2N + 1) * wigner3j(J, S, N, Ω, -Σ, -Λ)
#         η = (S = S, I = I, Λ = Λ, Σ = Σ, Ω = Ω, J = J, F = F, M = M)
#         a_state = HundsCaseB(; η...)

#         insert = true
#         for a_state_ in a_states
#             if a_state == a_state_
#                 insert = false
#             end
#         end
#         if insert  
#             push!(a_states, a_state)
#         end
#     end
