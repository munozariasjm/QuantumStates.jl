
module PhysicsTools
# Get the basis element with the maximum coefficient from the state
function get_max_coeff_basis_element(state)
    coeffs = abs.(state.coeffs) .^ 2
    max_idx = argmax(coeffs)
    return state.basis[max_idx]
end

# Extract quantum numbers from the basis element with the maximum coefficient
function extract_quantum_numbers_from_state(state)
    basis_elem = get_max_coeff_basis_element(state)
    quantum_numbers = Dict(:N => basis_elem.N, :J => basis_elem.J, :F => basis_elem.F, :M => basis_elem.M)
    return quantum_numbers
end

# Get the J quantum number from the basis element with the maximum coefficient
function get_J_quantum_number(state)
    basis_elem = get_max_coeff_basis_element(state)
    return basis_elem.J
end

# Search for a state in a list of states by specific quantum numbers
function find_state_by_quantum_numbers(states, N, J, F, M)
    for state in states
        qn = extract_quantum_numbers_from_state(state)
        if qn[:N] == N && qn[:J] == J && qn[:F] == F && qn[:M] == M
            return state
        end
    end
    error("State with quantum numbers N=$N, J=$J, F=$F, M=$M not found.")
end
export get_max_coeff_basis_element, extract_quantum_numbers_from_state, get_J_quantum_number, find_state_by_quantum_numbers
end
