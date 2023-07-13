using Parameters

Base.@kwdef struct HundsCaseA_LinearMolecule <: HundsCaseA
    E::Float64 = 0.0
    v_1::HalfInt
    v_2::HalfInt
    ℓ::HalfInt
    v_3::HalfInt
    Λ::HalfInt
    K::HalfInt # K = Λ + ℓ
    I::HalfInt
    S::HalfInt
    Σ::HalfInt
    J::HalfInt
    P::HalfInt # P = Λ + ℓ + Σ
    F::HalfInt
    M::HalfInt
    constraints = (
        K = Λ + ℓ,
        Σ = -S:S,
        # ℓ = unique(sign * (x + v2 % 2) for x ∈ 0:2:v_2 for sign ∈ (-1,1)),
        P = max(-J, Λ + ℓ + Σ):min(J, Λ + ℓ + Σ),
        F = abs(J - I):abs(J + I),
        M = -F:F
    )
end
export HundsCaseA_LinearMolecule

function unpack(state::HundsCaseA_LinearMolecule)
    (; v_1, v_2, ℓ, v_3, Λ, K, I, S, Σ, J, P, F, M) = state
    return v_1, v_2, ℓ, v_3, Λ, K, I, S, Σ, J, P, F, M
end
export unpack

function Rotation(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return (J * (J + 1) + S * (S + 1) + Λ^2 - 2 * P^2 + 2Λ * Σ) *
        δ(ℓ,ℓ′) * δ(Λ,Λ′) * δ(Σ,Σ′) * δ(J,J′) * δ(F,F′) * δ(M,M′)
end

function SpinOrbit(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return Λ * Σ * 
        δ(ℓ,ℓ′) * δ(Λ,Λ′) * δ(Σ,Σ′) * δ(J,J′) * δ(F,F′) * δ(M,M′)
end
export SpinOrbit

function SpinUncoupling(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    # Brown & Carrington, eq. (8.364)
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return -2 * (-1)^(J - P + S - Σ) * 
        sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
        sum(
            wigner3j_(J, 1, J, -P, q, P′) * 
            wigner3j_(S, 1, S, -Σ, q, Σ′)
            for q ∈ (-1,1)
        ) *
        δ(Λ,Λ′) * δ(ℓ,ℓ′) * δ(J,J′) * δ(F,F′) * δ(M,M′)
end
export SpinUncoupling

function ΛDoubling_q(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    """
    (Λ_+^2 J_-^2) term
    """
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return (
        (-1)^(J - P) *
        (1 / (2 * sqrt(6))) *
        sqrt( (2J - 1) * (2J) * (2J + 1) * (2J + 2) * (2J + 3) ) *
        sum(
            δ(Λ′, Λ + 2q) *
            wigner3j_(J, 2, J′, -P, -2q, P′)
            for q ∈ (-1,1)
        )
    ) * δ(ℓ,ℓ′) * δ(Σ,Σ′) * δ(F,F′) * δ(M,M′)
end
export ΛDoubling_q

function ΛDoubling_p2q(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    """
    (Λ_+^2 J_- S_-) term
    Brown and Carrington (eq. 9.66)
    """
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return (-1)^(J - P + S - Σ) * 
        sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
        sum(
            δ(Λ′, Λ + 2q) *
            wigner3j_(J, 1, J′, -P, -q, P′) *
            wigner3j_(S, 1, S′, -Σ, q, Σ′)
            for q ∈ (-1,1)
        ) *
        δ(ℓ,ℓ′) * δ(F,F′) * δ(M,M′)
    end
export ΛDoubling_p2q

function ℓDoubling(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    """
    H_ℓ = (1/2)q^v[(1/2)(G₊²J₋² + G₋²J₊²) - (G₊²J₋S₋ + G₋²J₊S₊)]
    """
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    term1 = (
        (-1)^(J - P) *
        (1 / (2 * sqrt(6))) *
        sqrt( (2J - 1) * (2J) * (2J + 1) * (2J + 2) * (2J + 3) ) *
        sum(
            δ(ℓ′, ℓ + 2q) *
            wigner3j_(J, 2, J′, -P, -2q, P′)
            for q ∈ (-1,1)
        )
    ) * δ(Σ,Σ′)
    # Add factors for G₊²???
    term2 = (
        (-1)^(J - P + S - Σ) * 
        sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
        sum(
            δ(ℓ′, ℓ + 2q) *
            wigner3j_(J, 1, J′, -P, -q, P′) *
            wigner3j_(S, 1, S′, -Σ, q, Σ′)
            for q ∈ (-1,1)
        )
    )
    return (term1 - term2) * δ(Λ,Λ′) * δ(F,F′) * δ(M,M′)
end
        
function Hyperfine_IL(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    # Brown and Carrington, eq. (8.372)
    # Orbital hyperfine interaction
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    if ~δ(Σ,Σ′) || ~δ(M,M′) || ~δ(F,F′) || ~δ(P,P′) || ~δ(Λ,Λ′)
        return 0.0
    else
        return (
            Λ * (-1)^(J′ + I + F + J - P) * 
            sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1)) *
            wigner6j_(I, J′, F, J, I, 1) * 
            wigner3j_(J, 1, J′, -P, 0, P′) * 
            δ(Λ, Λ′) * δ(F, F′) * δ(M, M′)
        )
    end
end
export Hyperfine_IL

function Hyperfine_IF(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    # Fermi contact interaction
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return (-1)^(I + J′ + F + S - Σ + J - P) * 
        sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1)) *
        wigner6j_(I, J′, F, J, I, 1) *
        sum(
            wigner3j_(J, 1, J′, -P, q, P′) *
            wigner3j_(S, 1, S, -Σ, q, Σ′)
            for q in -1:1
        ) *
        δ(F, F′) * δ(M, M′)
end
export Hyperfine_IF
    
function Hyperfine_Dipolar_c(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return sqrt(30) * (1/3) * (-1)^(I + J′ + F) * (-1)^(J - P) * (-1)^(S - Σ) *
        wigner6j_(I, J′, F, J, I, 1) *
        sqrt( I * (I + 1) * (2I + 1) ) *
        sqrt( (2J + 1) * (2J′ + 1) ) *
        sqrt( S * (S + 1) * (2S + 1) ) *
        sum(
            (-1)^q * 
            wigner3j_(J, 1, J′, -P, q, P) *
            sum(
                wigner3j_(1, 2, 1, q′, 0, -q) *
                wigner3j_(S, 1, S, -Σ, q′, Σ′)
                for q′ ∈ -1:1
            ) for q ∈ -1:1
        ) * 
        δ(F, F′) * δ(M, M′)
end
export Hyperfine_Dipolar_c

function Hyperfine_Dipolar_d(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return sqrt(30) * (1/2) * sqrt(3/2) * (2/3) * (-1)^(I + J′ + F) * (-1)^(J - P) * (-1)^(S - Σ) *
        wigner6j_(I, J′, F, J, I, 1) *
        sqrt( I * (I + 1) * (2I + 1) ) *
        sqrt( (2J + 1) * (2J′ + 1) ) *
        sqrt( S * (S + 1) * (2S + 1) ) *
        sum(
            (-1)^q * 
            wigner3j_(J, 1, J′, -P, q, P) *
            sum(
                (
                    wigner3j_(1, 2, 1, q′, 2, -q) +
                    wigner3j_(1, 2, 1, q′, -2, -q)
                ) *
                wigner3j_(S, 1, S, -Σ, q′, Σ′)
                for q′ ∈ -1:1
            ) for q ∈ -1:1
        ) * 
        δ(F, F′) * δ(M, M′)
end
export Hyperfine_Dipolar_d

function RennerTeller(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    """
    See Brown (1975).
    """
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return (
        sum(((v_2 + 1)^2 - K^2)^(1/2) * δ(Λ,Λ′+2q) * δ(ℓ,ℓ′-2q) for q ∈ (-1,1))
    ) * 
    δ(Σ,Σ′) * δ(J,J′) * δ(F,F′) * δ(M,M′)
end
export RennerTeller

function TDM(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule, p::Integer)
    v_1,  v_2,  ℓ,  v_3,  Λ,  K,  I,  S,  Σ,  J,  P,  F,  M  = unpack(state)
    v_1′, v_2′, ℓ′, v_3′, Λ′, K′, I′, S′, Σ′, J′, P′, F′, M′ = unpack(state′)
    return (-1)^(F′ - M′) * wigner3j(F′, 1, F, -M′, p, M) * (-1)^(F + J′ + I + 1) * sqrt((2F′ + 1) * (2F + 1)) *
        δ(Σ, Σ′) * (-1)^(J′ - P′) * wigner6j(J, F, I, F′, J′, 1) * sqrt((2J′ + 1) * (2J + 1)) *
        sum(
            wigner3j(J′, 1, J, -P′, q, P) for q ∈ -1:1
        )
end
TDM(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule) = sum(TDM(state, state′, p) for p ∈ -1:1)
export TDM

# function Zeeman_L(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
#     """
#     Electron orbit Zeeman
#     """
#     I,  S,  Λ,  J,  P,  Σ,  F,  M  = unpack(state)
#     I′, S′, Λ′, J′, P′, Σ′, F′, M′ = unpack(state′)

#     if ~δ(P,P′) || ~δ(Σ,Σ′) || ~δ(Λ,Λ′) || ~δ(M,M′) || ~δ(S,S′) || ~δ(I,I′)
#         return 0.0
#     else
#         return (
#             Λ * (-1)^(F-M+F′+J+I+1+J-P) * wigner6j(J,F,I,F′,J′,1) * 
#             wigner3j(F,1,F′,-M,0,M) * sqrt((2F+1)*(2F′+1)*(2J+1)*(2J′+1)) * 
#             wigner3j(J,1,J′,-P,0,P)
#         )
#     end
# end

# function Zeeman_S(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
#     """
#     Electron spin Zeeman
#     """
#     I,  S,  Λ,  J,  P,  Σ,  F,  M  = unpack(state)
#     I′, S′, Λ′, J′, P′, Σ′, F′, M′ = unpack(state′)

#     if ~δ(S,S′) || ~δ(I,I′) || ~δ(M,M′) || ~δ(Λ,Λ′)
#         return 0.0
#     else
#         return (
#                 sum(
#                     (-1)^(F-M+J+I+F′+1+J-P+S-Σ) * wigner6j(J,F,I,F′,J′,1) * 
#                     wigner3j(F,1,F′,-M,0,M′) * sqrt((2F+1)*(2F′+1)) * 
#                     wigner3j(J,1,J′,-P,q,P′) * sqrt((2J+1)*(2J′+1)) *
#                     wigner3j(S,1,S,-Σ,q,Σ′) * sqrt(S*(S+1)*(2S+1))
#                     for q ∈ -1:1
#                 )
#             )
#     end
# end

# function Zeeman_gl(state::HundsCaseA_LinearMolecule, state′::HundsCaseA_LinearMolecule)
#     """
#     From CaH optical Zeeman paper, Eq. 4
#     """
#     I,  S,  Λ,  J,  P,  Σ,  F,  M  = unpack(state)
#     I′, S′, Λ′, J′, P′, Σ′, F′, M′ = unpack(state′)

#     if ~δ(Λ,Λ′) || ~δ(S,S′) || ~δ(I,I′)
#         return 0.0
#     else
#         return (
#             (-1)^(F-M) * wigner3j_(F,1,F′,-M,0,M′) * (-1)^(F′+J+I+1) * sqrt((2F+1)*(2F′+1)) *
#             wigner6j_(J′,F′,I,F,J,1) * 
#             sum(
#                 (-1)^(J-P) * wigner3j_(J,1,J′,-P,q,P′) * sqrt((2J+1)*(2J′+1)) *
#                 (-1)^(S-Σ) * wigner3j_(S,1,S,-Σ,q,Σ′) * sqrt((2S+1)*(2S′+1))
#                 for q ∈ (-1,1)
#             )
#         )
#     end
# end