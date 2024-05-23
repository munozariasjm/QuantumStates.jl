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

function overlap(state::HundsCaseB, state′::HundsCaseA)
    # Eq. (6.149) in Brown & Carrington, note the equation has an error
    S,  I,  Λ,  N, J, F,  M      = state.S, state.I, state.Λ, state.N, state.J, state.F, state.M
    S′, I′, Λ′, Σ, Ω, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    return (-1)^(N - S + Ω) * sqrt(2N + 1) * wigner3j(J, S, N, Ω, -Σ, -Λ) * δ(J, J′) * δ(F, F′) * δ(M, M′) * δ(Λ, Λ′)
end
overlap(state::HundsCaseA, state′::HundsCaseB) = overlap(state′, state)

function overlap(state::HundsCaseB, state′::HundsCaseA_LinearMolecule)
    # Eq. (6.149) in Brown & Carrington, note the equation has an error
    S,  I,  ℓ,  N, J, F,  M      = state.S, state.I, state.Λ, state.N, state.J, state.F, state.M
    S′, I′, ℓ′, Σ, P, J′, F′, M′ = state′.S, state′.I, state′.ℓ, state′.Σ, state′.P, state′.J, state′.F, state′.M
    return (-1)^(-J + P + 2S) * sqrt(2N + 1) * wigner3j(J, N, S, P, -ℓ, -Σ) * δ(J, J′) * δ(F, F′) * δ(M, M′)
end
overlap(state::HundsCaseA_LinearMolecule, state′::HundsCaseB) = overlap(state′, state)

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

function Parity(state::State{<:HundsCaseB})
    # Note: For an actual parity state, only two basis states should have nonzero coefficients
    basis_state_idxs = findall(c -> norm(c) > 0.45, state.coeffs)
    basis_states = state.basis[basis_state_idxs]
    
    @unpack N, Λ = basis_states[1]
    sign1 = sign(state.coeffs[basis_state_idxs[1]])
    phase = (-1)^Λ * sign1
    @unpack N, Λ = basis_states[2]
    sign2 = sign(state.coeffs[basis_state_idxs[2]])
    
    return -(-1)^(N-abs(Λ)) * (-1)^(sign1 == sign2)
end
export Parity

# Added 5/2/24
function overlap(state::HundsCaseB_LinearMolecule, state′::HundsCaseA_LinearMolecule)
    # See Hirota eq. (2.3.3)
    S,  I,  ℓ, Λ, K, N, J, F, M = state.S, state.I, state.ℓ, state.Λ, state.K, state.N, state.J, state.F, state.M
    S′, I′, ℓ′, Λ′, Σ, P, K′, J′, F′, M′ = state′.S, state′.I, state′.ℓ, state′.Λ, state′.Σ, state′.P, state′.K, state′.J, state′.F, state′.M
    return (-1)^(-J + P + 2S) * sqrt(2N + 1) * wigner3j(J, N, S, P, -K, -Σ) * δ(J,J′) * δ(F,F′) * δ(M,M′) * δ(ℓ,ℓ′) * δ(Λ,Λ′)
end
overlap(state::HundsCaseA_LinearMolecule, state′::HundsCaseB_LinearMolecule) = overlap(state′, state)