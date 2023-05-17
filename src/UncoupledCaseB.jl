using Parameters

# abstract type HundsCaseB_Decoupled <: BasisState end
# export HundsCaseB_Decoupled

@composite Base.@kwdef struct HundsCaseB_Decoupled <: BasisState
    E::Float64 = 0.0
    Λ::HalfInt
    N::HalfInt
    M_N::HalfInt
    S::HalfInt
    M_S::HalfInt
    I::HalfInt
    M_I::HalfInt
    constraints = (
        N = abs(Λ):∞,
        M_N = -N:N,
        M_S = -S:S,
        M_I = -I:I
    )
end
export HundsCaseB_Decoupled

function overlap(state::HundsCaseB, state′::HundsCaseB_Decoupled)
    @unpack S, I, Λ, N, J, F, M = state
    Λ′ = Λ
    N′ = N
    M_F = M
    @unpack Λ, N, M_N, S, M_S, I, M_I = state′

    return δ(Λ, Λ′) * δ(N, N′) * sum(
          (-1)^(J - I + M_F) * sqrt(2F + 1) * wigner3j(J, I, F, M_J, M_I, -M_F)
        * (-1)^(N - S + M_J) * sqrt(2J + 1) * wigner3j(N, S, J, M_N, M_S, -M_J)
        for M_J ∈ -J:J
    )
end

overlap(state′::HundsCaseB_Decoupled, state::HundsCaseB) = overlap(state, state′)
export overlap