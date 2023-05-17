struct TensorProductState{T1<:BasisState,T2<:BasisState} <: BasisState
    basis_state1::T1
    basis_state2::T2
end
export TensorProductState

overlap(state::TensorProductState{T1,T2}, state′::TensorProductState{T3,T4}) where {T1,T2,T3,T4} = 
    overlap(state.basis_state1, state′.basis_state1) * overlap(state.basis_state2, state′.basis_state2)