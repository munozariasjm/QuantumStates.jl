Base.@kwdef struct AngularMomentumState_withSpinRotation_Uncoupled <: BasisState
    E::Float64 = 0.0
    N::HalfInt
    M_N::HalfInt
    S::HalfInt
    M_S::HalfInt
    constraints = (
        M_N = -N:N,
        M_S = -S:S
    )
end
export AngularMomentumState_withSpinRotation_Uncoupled

function unpack(state::AngularMomentumState_withSpinRotation_Uncoupled)
    return (state.N, state.M_N, state.S, state.M_S)
end
export unpack

function overlap(state::AngularMomentumState_withSpinRotation, state′::AngularMomentumState_withSpinRotation_Uncoupled)
    J, N′, S, M_J  = unpack(state)
    N, M_N, S, M_S = unpack(state′)
    # return δ(N, N′) * sum(
    #     (-1)^(N - S + M_J) * sqrt(2J + 1) * wigner3j(N, S, J, M_N, M_S, -M_J)
    #     for M_N ∈ -N:N, M_S ∈ -S:S
    # )
    return δ(N, N′) * (-1)^(N - S + M_J) * sqrt(2J + 1) * wigner3j(N, S, J, M_N, M_S, -M_J)
end
overlap(state′::AngularMomentumState_withSpinRotation_Uncoupled, state::AngularMomentumState_withSpinRotation) = overlap(state, state′)
export overlap

# function overlap(state::AngularMomentumState_withSpinRotation, state′::AngularMomentumState_withSpinRotation_Uncoupled)
#         J, N′, S, M_J  = unpack(state)
#         N, M_N, S, M_S = unpack(state′)
#         # return δ(N, N′) * sum(
#         #     (-1)^(N - S + M_J) * sqrt(2J + 1) * wigner3j(N, S, J, M_N, M_S, -M_J)
#         #     for M_N ∈ -N:N, M_S ∈ -S:S
#         # )
#         return δ(N, N′) * sum(
#             (-1)^(N - S + Ω) * sqrt(2N + 1) * wigner3j(N, S, J, M_N, M_S, -Ω)
#             for Ω ∈ -J:J)
#     end
#     overlap(state′::AngularMomentumState_withSpinRotation_Uncoupled, state::AngularMomentumState_withSpinRotation) = overlap(state, state′)
#     export overlap