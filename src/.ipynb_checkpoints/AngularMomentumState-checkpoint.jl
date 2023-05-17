using Parameters
using InfiniteArrays
using CompositeStructs
using HalfIntegers

Base.@kwdef mutable struct AngularMomentumState <: BasisState
    E::Float64 = 0.0
    F::HalfInt
    M::HalfInt
    constraints = (
        M = -F:F,
    )
end
export AngularMomentumState

function unpack(state::AngularMomentumState)
    return (state.F, state.M)
end
export unpack

function TDM(state::AngularMomentumState, state′::AngularMomentumState, p::Int64)
    F,  M  = unpack(state)
    F′, M′ = unpack(state′)
    return (
        (-1)^p * (-1)^(F - M) * wigner3j_(F, 1, F′, M, p, -M′) * sqrt(2F + 1)
    )
end

function TDM_magnetic(state::AngularMomentumState, state′::AngularMomentumState, p::Int64)
    # Assumes magnetic moment aligned along z-axis of molecule-fixed axis
    F,  M  = unpack(state)
    F′, M′ = unpack(state′)
    return -(
        (-1)^p * (-1)^(F - M) * sqrt(F * (F + 1) * (2F + 1)) * wigner3j(F, 1, F′, M, p, -M′)
    )
end
# TDM_magnetic(state::State, state′::State, p::Int64) = extend_operator(TDM_magnetic, state, state′, p)

function H(state::AngularMomentumState, state′::AngularMomentumState)
    δ(state, state′) * state.E
end
export H