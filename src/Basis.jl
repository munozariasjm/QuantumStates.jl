abstract type BasisState end
export BasisState

"""
    Two basis states are equal if all quantum numbers are the same.
"""
# function ==(state::BasisState, state′::BasisState)
#     state_type = typeof(state)
#     for field in fieldnames(state_type)
#         if getfield(state, field) != getfield(state′, field)
#             return false
#         end
#     end
#     return true
# end
# export ==

function enumerate_states(state_type, QN_bounds)
    states = state_type[]
    QNs = [QN for QN ∈ fieldnames(state_type) if QN ∉ (:E, :constraints)]
    
    # Define a state with all QN = 0 to get the constraints for the QNs
    η = NamedTuple([QN => 0 for QN in QNs])
    
    enumerate_states(η, states, state_type, QNs, QN_bounds, 1)

    return states
end

function enumerate_states(η, states, state_type, QNs, QN_bounds, idx, max_states=1000)

    if length(states) > max_states
        return false
    end

    if idx == length(QNs) + 1
        new_state = state_type(; η...)
        push!(states, new_state)
        return nothing
    end
    
    # Check if quantum number has been given bounds; else apply constraints
    iterated_QN = QNs[idx]
    QN_constraints = state_type(; η...).constraints
    if iterated_QN ∈ keys(QN_constraints)
        QN_constraint = QN_constraints[iterated_QN]
        QN_constraint_bounds = eval(QN_constraint)
    end
    
    if iterated_QN ∈ keys(QN_bounds)
        bounds = QN_bounds[iterated_QN]
        for i ∈ bounds
            if iterated_QN ∉ keys(QN_constraints) || (i ∈ QN_constraint_bounds)
                η′ = (; η..., iterated_QN => i)
                enumerate_states(η′, states, state_type, QNs, QN_bounds, idx + 1)
            end
        end
    elseif iterated_QN ∈ keys(QN_constraints)
        for i ∈ QN_constraint_bounds
            η′ = (; η..., iterated_QN => i)
            enumerate_states(η′, states, state_type, QNs, QN_bounds, idx + 1)
        end
    else
        enumerate_states(η, states, state_type, QNs, QN_bounds, idx + 1)
    end
    
    return nothing
end
export enumerate_states

# Need to redefine so that this is type-stable...
# function enumerate_states(state_type, QN_bounds1, QN_bounds2)
#     basis_states1 = enumerate_states(state_type.types[1], QN_bounds1)
#     basis_states2 = enumerate_states(state_type.types[2], QN_bounds2)
#     basis = state_type[]
#     for basis_state1 ∈ basis_states1
#         for basis_state2 ∈ basis_states2
#             push!(basis, state_type(basis_state1, basis_state2))
#         end
#     end
#     return basis
# end