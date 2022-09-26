using Parameters
using InfiniteArrays
using CompositeStructs
using HalfIntegers

abstract type Atom <: BasisState end
export Atom

@composite Base.@kwdef struct Atom
    E::Float64 = 0.0
    L::HalfInt
    S::HalfInt
    J::HalfInt
    I::HalfInt
    F::HalfInt
    M::HalfInt
    constraints = (
        J = abs(L - S):abs(L + S),
        F = abs(J - I):abs(J + I),
        M = -F:F
    )
end
export Atom

function unpack(state::Atom)
    return (state.L, state.S, state.J, state.I, state.F, state.M)
end
export unpack

function TDM(state::Atom, state′::Atom, p::Int64)
    L,  S,  J,  I,  F,  M = unpack(state)
    L′, S′, J′, I′, F′, M′ = unpack(state′)
    return (
          (-1)^p * (-1)^(F - M) * wigner3j(F, 1, F′, -M, p, M′)
        * (-1)^(J + I + F′ + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * wigner6j(J, F, I, F′, J′, 1)
        * (-1)^(L + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) ) * wigner6j(L, J, S, J′, L′, 1)
        * sqrt(2L + 1)
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
