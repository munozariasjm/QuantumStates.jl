using Parameters
using InfiniteArrays
using CompositeStructs
using HalfIntegers

@composite Base.@kwdef struct AsymmetricTopMolecule <: HundsCaseB
    E::Float64 = 0.0
    v1::HalfInt
    v2::HalfInt
    v3::HalfInt
    S::HalfInt
    I::HalfInt
    N::HalfInt
    K_A::HalfInt
    K_C::HalfInt
    J::HalfInt
    F::HalfInt
    M::HalfInt
    constraints = (
        K_A = 0:N,
        K_C = (N - K_A):N,
        J = abs(N - S):abs(N + S),
        F = abs(J - I):abs(J + I),
        M = -F:F
    )
end
export AsymmetricTopMolecule

function Rotation_A(state::HundsCaseB, state′::HundsCaseB)
    S, I, Λ, N, J, F, M = unpack(state)
    S′, I′, Λ′, N′, J′, F′, M′ = unpack(state′)
    if ~δ(Λ, Λ′) || ~δ(N, N′) || ~δ(J, J′) || ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else
        return N * (N + 1)
    end
end
export Rotation