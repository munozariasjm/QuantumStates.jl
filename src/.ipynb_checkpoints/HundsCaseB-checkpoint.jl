using Parameters
using InfiniteArrays
using CompositeStructs
using HalfIntegers

# Define the spherical tensor T^k_q(ϵ), here for linear and symmetric top molecules
const T = [
    0.0 0.0 0.0
    -2/√3 0.0 -2/√6
    0.0 0.0 0.0
    ]
export T

abstract type HundsCaseB <: BasisState end
export HundsCaseB

@composite Base.@kwdef struct HundsCaseB_Rot <: HundsCaseB
    E::Float64 = 0.0
    S::HalfInt 
    I::HalfInt
    Λ::HalfInt
    N::HalfInt
    J::HalfInt
    F::HalfInt
    M::HalfInt
    constraints = (
        N = abs(Λ):∞,
        J = abs(N - S):abs(N + S),
        F = abs(J - I):abs(J + I),
        M = -F:F
    )
end
export HundsCaseB_Rot

function unpack(state::HundsCaseB)
    return (state.S, state.I, state.Λ, state.N, state.J, state.F, state.M)
end
export unpack

function Rotation(state::HundsCaseB, state′::HundsCaseB)
    S, I, Λ, N, J, F, M = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(Λ, Λ′) || ~δ(N, N′) || ~δ(J, J′) || ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else
        return N * (N + 1) - Λ^2
    end
end
export Rotation

function RotationDistortion(state::HundsCaseB, state′::HundsCaseB)
    S, I, Λ, N, J, F, M = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(Λ, Λ′) || ~δ(N, N′) || ~δ(J, J′) || ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else
        return - (N * (N + 1) - Λ^2)^2
    end
end
export RotationDistortion
    
# Spin-rotation for zero internuclear axis angular momentum, i.e., Λ = 0
function SpinRotation_Λ0(state::HundsCaseB, state′::HundsCaseB)
    S, I, Λ, N, J, F, M = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(Λ, Λ′) || ~δ(N, N′) || ~δ(J, J′) || ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else 
        return (
                (-1)^(N + S + J) 
                * sqrt(S * (S + 1) * (2S + 1) * N * (N + 1) * (2N + 1))
                * wigner6j(S, N, J, N, S, 1)
            )
    end
end
export SpinRotation_Λ0

# Spin-rotation for Λ != 0, reduces to above matrix element for Λ = 0
function SpinRotation(state::HundsCaseB, state′::HundsCaseB)
    # Hirota, eq. (2.3.35)
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(J, J′) || ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else
        return (
                (1/2) * (-1)^(J′ + S + N) * (-1)^(N - Λ) 
                * sqrt( S * (S + 1) * (2S + 1) * (2N + 1) * (2N′ + 1) ) 
                * wigner6j(N′, S ,J, S, N, 1)
                * sum( sqrt(2k + 1)
                * (
                    (-1)^k * 
                    sqrt( N′ * (N′ + 1) * (2N′ + 1) ) * wigner6j(1, 1, k, N, N′, N′) 
                    + 
                    sqrt( N * (N + 1) * (2N + 1) ) * wigner6j(1, 1, k, N′, N, N)
                    )
                    * wigner3j(N, k, N′, -Λ, q, Λ′) * T[q+2, k+1]
                for k in 0:2, q in -1:1
            )
        )
    end
end
export SpinRotation

function Hyperfine_IS(state::HundsCaseB, state′::HundsCaseB)
    # Fermi-contact interaction
    # Hirota, pg. 39
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(Λ, Λ′) || ~δ(N, N′) || ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else
        return (
                (-1)^(N′ + S + J) * (-1)^(J′ + I + F′ + 1)
                * sqrt( (2J′ + 1) * (2J + 1) * S * (S + 1) * (2S + 1) * I * (I + 1) * (2I + 1) )
                * wigner6j(I, J, F′, J′, I, 1)
                * wigner6j(S, J, N′, J′, S, 1)
            )
    end
end
export Hyperfine_IS

function Hyperfine_Dipolar(state::HundsCaseB, state′::HundsCaseB)
    # Dipolar interaction term, from c(Iz ⋅ Sz)
    # Hirota, pg. 39
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else
        return (
                sqrt(30) * (-1)^(N - Λ) * (-1)^(J′ + I + F + 1) *
                wigner6j(I, J, F′, J′, I, 1) * 
                wigner9j(N, N′, 2, S, S, 1, J, J′, 1) * 
                wigner3j(N, 2, N′, -Λ, 0, Λ′) *
                sqrt( S * (S + 1) * (2S + 1) * I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * (2N + 1) * (2N′ + 1) )
            )
    end
end
export Hyperfine_Dipolar

function ℓDoubling(state::HundsCaseB, state′::HundsCaseB)
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(N, N′) || ~δ(J, J′) || ~δ(F, F′) || ~δ(M, M′) || ~δ(abs(Λ′ - Λ), 2)
        return 0.0
    else
        return (
                (-1)^(N - Λ) *
                (1 / (2 * sqrt(6))) *
                sqrt( (2N - 1) * (2N) * (2N + 1) * (2N + 2) * (2N + 3) ) *
                sum(
                    wigner3j(N, 2, N′, -Λ, 2q, Λ′)
                    for q ∈ (-1,1)
                )                   
            )
    end
end
export ℓDoubling

# function Hyperfine_IK(state::HundsCaseB, state′::HundsCaseB)
#     S, I, N, Λ, J, F, M = unpack(state)   
#     S′, I′, N′, Λ′, J′, F′, M′ = unpack(state′)
#     return (-1)^(F + I + N′ + S + 2J + 1 + N - Λ) * 
#         sqrt( I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * (2N + 1) * (2N′ + 1) ) *
#         wigner6j(I, J, F, J′, I, 1) * 
#         wigner6j(N, J, S, J′, N′, 1) * 
#         wigner3j(N′, 1, N, -Λ, 0, Λ) *
#         δ(Λ, Λ′) * δ(F, F′) * δ(M, M′)
# end
# export Hyperfine_IK

# function Hyperfine_SK(state::HundsCaseB, state′::HundsCaseB)
#     S, I, N, Λ, J, F, M = unpack(state)   
#     S′, I′, N′, Λ′, J′, F′, M′ = unpack(state′)
#     return (-1)^(2N + J + S - Λ) * sqrt(S * (S + 1) * (2S + 1) * (2N + 1) * (2N′ + 1)) *
#         wigner6j(S, N, J, N′, S, 1) * 
#         wigner3j(N′, 1, N, -Λ, 0, Λ) *
#         δ(Λ, Λ′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
# end
# export Hyperfine_SK

function Stark(state::HundsCaseB, state′::HundsCaseB)
    # Hirota, equation (2.5.35)
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    return (
            -(-1)^(F - M) * wigner3j(F, 1, F′, -M, 0, M′) 
            * (-1)^(J + I + F′ + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner6j(J, F, I, F′, J′, 1)
            * (-1)^(N + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) ) * wigner6j(N, J, S, J′, N′, 1)
            * (-1)^(N - Λ) * sqrt( (2N + 1) * (2N′ + 1) ) * wigner3j(N, 1, N′, -Λ, 0, Λ′)
    )
end
export Stark

function Zeeman(state::HundsCaseB, state′::HundsCaseB, p::Int64)
    # Hirota, equation (2.5.16) and (2.5.19)
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(Λ, Λ′) || ~δ(N, N′)
        return 0.0
    else
        return (
                  (-1)^p * (-1)^(F′ - M′) * wigner3j(F′, 1, F, -M′, -p, M)
                * (-1)^(J′ + I + F + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner6j(J′, F′, I, F, J, 1)
                * (-1)^(N + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1) ) * wigner6j(S, J′, N, J, S, 1)
        )
    end
end
export Zeeman

# function Zeeman(state::HundsCaseB, state′::HundsCaseB, B::Vector{Float64})
#     # Hirota, equation (2.5.16) and (2.5.19)
#     S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
#     S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
#     if ~δ(Λ, Λ′) || ~δ(N, N′)
#         return 0.0
#     else
#         return (
#                   (-1)^(J′ + I + F + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner6j(J′, F′, I, F, J, 1)
#                 * (-1)^(N + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1) ) * wigner6j(S, J′, N, J, S, 1)
#             ) * 
#             sum(
#                 B[p+2] * (-1)^p * (-1)^(F′ - M′) * wigner3j(F′, 1, F, -M′, -p, M) for p ∈ -1:1
#             )
#     end
# end
# export Zeeman

function Σ(state::HundsCaseB)
    @unpack Λ, N, S, J = state
    val = zero(Float64)
    for Σ ∈ -S:S
        Ω = Λ + Σ
        val += Σ * (2N + 1) * wigner3j(J, S, N, Ω, -Σ, -Λ)^2
    end
    return val
end
# function Σ(state::State)
#     val = zero(Float64)
#     for i ∈ eachindex(state.basis), j ∈ eachindex(state.basis)
#         val += conj(state.coeffs[i]) * Σ(state.basis[j]) * state.coeffs[j]
#     end
#     return val
# end
Σ(state::State) = sum(Σ(state.basis[i]) * state.coeffs[i] * conj(state.coeffs[i]) for i ∈ eachindex(state.basis))
export Σ

function TDM_magnetic(state::HundsCaseB, state′::HundsCaseB, p::Int64)
    # Assumes magnetic moment aligned along z-axis of molecule-fixed axis
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(Λ, Λ′) || ~δ(N, N′)
        return 0.0
    else
        return (
                (-1)^p * (-1)^(F′ - M′) * wigner3j(F′, 1, F, -M′, -p, M)
                * (-1)^(J′ + I + F + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner6j(J′, F′, I, F, J, 1)
                * (-1)^(N + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1) ) * wigner6j(S, J′, N, J, S, 1)
            )
    end
end
TDM_magnetic(state::State, state′::State, p::Int64) = extend_operator(TDM_magnetic, state, state′, p)
export TDM_magnetic

function TDM(state::HundsCaseB, state′::HundsCaseB, p::Int64)
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    return (
          (-1)^p * (-1)^(F - M) * wigner3j(F, 1, F′, -M, p, M′)
        * (-1)^(J + I + F′ + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner6j(J, F, I, F′, J′, 1)
        * (-1)^(N + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) ) * wigner6j(N, J, S, J′, N′, 1)
        * (-1)^(N - Λ) * sqrt( (2N + 1) * (2N′ + 1) ) * sum(wigner3j(N, 1, N′, -Λ, q, Λ′) for q ∈ -1:1)
    )
end
TDM(state, state′) = sum(TDM(state, state′, p) for p ∈ -1:1)

function 𝒫(K,P,ϵ)
    val = 0.0
    ϵm1, ϵ0, ϵp1 = ϵ
    if P == 0
        if K == 0
            val += 1.0
        elseif K == 1
            val += ϵp1 * conj(ϵp1) - ϵm1 * conj(ϵm1)
        elseif K == 2
            val += -(1/2) * (1 - 3 * ϵ0 * conj(ϵ0))
        end
    elseif P == +1
        if K == 1
            val += - (ϵ0 * conj(ϵm1) + conj(ϵ0) * ϵp1)
        elseif K == 2
            val += sqrt(3/2) * (-ϵ0 * conj(ϵm1) + conj(ϵ0) * ϵp1)
        end
    elseif P == -1
        if K == 1
            val += (ϵ0 * conj(ϵp1) + conj(ϵ0) * ϵm1)
        elseif K == 2
            val += sqrt(3/2) * (-ϵ0 * conj(ϵp1) + conj(ϵ0) * ϵm1)
        end
    elseif P == +2
        if K == 2
            val += -sqrt(3/2) * conj(ϵm1) * ϵp1
        end
    elseif P == -2
        if K == 2
            val += -sqrt(3/2) * conj(ϵp1) * ϵm1
        end
    end
    return val
end             
export 𝒫
                                                            
function polarizability(state::HundsCaseB, state′::HundsCaseB, α, ϵ)
    S,  I,  Λ,  N,  J,  F,  M  = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    val = 0.0
    for K in 0:2
        for P in -K:K
#             val += (
#                 (-1)^P
#                 * (-1)^(F′ - M′) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner3j(F′, K, F, -M′, P, M)
#                 * (-1)^(F + J′ + I) * wigner6j(J′, F′, I, F, J, K)
#                 * ((-1)^(N′ + N) + 1) * sqrt( (2N + 1) * (2N′ + 1) )
#                 * wigner3j(J,  S,  N,  -half(1) + Λ,  half(1), -Λ) 
#                 * wigner3j(J′, S′, N′, -half(1) + Λ′, half(1), -Λ′)
#                 * (-1)^(J + S) * sqrt( (2J + 1) * (2J′ + 1) ) 
# #                 * wigner3j(J′, K, J, -half(1) - Λ′, (Λ′ - Λ), half(1) + Λ)
#                 * wigner3j(J′, K, J, -half(1) - Λ′, 0, half(1) + Λ)
#                 * α[K+1] * 𝒫(K, -P, ϵ)
#                 * (-1)^(Λ + Λ′)
#             )
            val += (
#                 * (-1)^(F′ - M′) * wigner3j(F′, K, F, -M′, P, M)
#                 * (-1)^(F + J′ + K + I) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner6j(J′, F′, I, F, J, K)
#                 * (-1)^(N′ - Λ′) * sqrt( (2N + 1) * (2N′ + 1) ) * wigner3j(N′, K, N, -Λ′, 0, Λ)
#                 * (-1)^(J + N′ + K + S) * sqrt( (2J + 1) * (2J′ + 1) ) * wigner6j(N′, J′, S, J, N, K)
                -(-1)^P
                * (-1)^(F - M)
                * wigner3j(F, K, F′, -M, P, M′)
                * (-1)^(J + I + F′ + K) * sqrt( (2F + 1) * (2F′ + 1) )
                * wigner6j(J, F, I, F′, J′, K)
                * (-1)^(N + S + J′ + K) * sqrt( (2J + 1) * (2J′ + 1) )
                * wigner6j(N, J, S, J′, N′, K)
                * (-1)^(N - Λ) * sqrt( (2N + 1) * (2N′ + 1) )
                * wigner3j(N, K, N′, -Λ, 0, Λ′)
                * α[K+1] * 𝒫(K, -P, ϵ)
                # * δ(Λ, Λ′)
            ) 
        end
    end
    return val
end
export polarizability
