using Parameters
using InfiniteArrays
using CompositeStructs
using HalfIntegers

@composite Base.@kwdef struct AsymmetricTopMolecule <: HundsCaseB
    E::Float64 = 0.0
    v_1::HalfInt
    v_2::HalfInt
    v_3::HalfInt
    S::HalfInt
    I::HalfInt
    N::HalfInt
    K::HalfInt
    J::HalfInt
    F::HalfInt
    M::HalfInt
    constraints = (
        K = -N:N,
        J = abs(N - S):abs(N + S),
        F = abs(J - I):abs(J + I),
        M = -F:F
    )
end
export AsymmetricTopMolecule

# Define the tensor for the rotational constants
# Based on Comaker et al. (1973)
const T_Rotational = [
    0 0 -1/√3 0 0
    0 0 0 0 0
    1/√6 0 0 0 1/2
    ]
const T_A = T_Rotational .* [
    0 0 1 0 0
    0 0 0 0 0
    0 0 2 0 0
    ]
const T_B = T_Rotational .* [
    0 0 1 0 0 
    0 0 0 0 0
    0 -1 0 1 0
    ]
const T_C = T_Rotational .* [
    0 0 1 0 0 
    0 0 0 0 0
    0 -1 0 -1 0
    ];

function Rotation(state::AsymmetricTopMolecule, state′::AsymmetricTopMolecule, T)
    N,  J,  F,  M,  K  = state.N,  state.J,  state.F,  state.M,  state.K
    N′, J′, F′, M′, K′ = state′.N, state′.J, state′.F, state′.M, state′.K
    if ~δ(N, N′) || ~δ(J, J′) || ~δ(F, F′) || ~δ(M, M′)
        return 0.0
    else
        return (
            (-1)^(N′ - K) * N * (N + 1) * (2N + 1) *
            sum(
                T[k+1,q+3] * wigner3j(N, k, N, -K, q, K′)
                * (-1)^k * (2k + 1)^(1/2) * wigner6j(N, N, 1, k, 1, N)
                for k ∈ 0:2, q ∈ -2:2
            )
        )
    end
end
Rotation_A(state, state′) = Rotation(state, state′, T_A)
Rotation_B(state, state′) = Rotation(state, state′, T_B)
Rotation_C(state, state′) = Rotation(state, state′, T_C)
export Rotation, Rotation_A, Rotation_B, Rotation_C
