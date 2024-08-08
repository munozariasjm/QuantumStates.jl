# -*- coding: utf-8 -*-
using Pkg
Pkg.activate("/Users/jose/Documents/Works/MIT/RaX/Simu/Molecule-Sims/")
using
    Revise,
    QuantumStates
;
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

function define_YbF_states(B)

    # #################################################################################################
    # Ground state Hamiltonian
    # #################################################################################################
    QN_bounds = (
        S=1 / 2,           # Spin quantum number
        I=1 / 2,           # Nuclear spin quantum number
        Λ=0,             # Projection of the electronic orbital angular momentum
        N=0:3            # Rotational quantum number range
    )
    X_state_basis = enumerate_states(HundsCaseB_LinearMolecule, QN_bounds)
    X_state_operator = :(
        BX * Rotation +                     # Rotational energy term
        DX * RotationDistortion +           # Rotational distortion energy term
        γX * SpinRotation +                 # Spin-rotation interaction term
        bFX * Hyperfine_IS +                # Hyperfine interaction (Fermi contact term)
        cX * (Hyperfine_Dipolar / 3)        # Hyperfine interaction (dipolar term)
    )

    X_state_parameters = QuantumStates.@params begin
        BX = 7233.8271e6         # Rotational constant in Hz
        DX = 0.0                  # Rotational distortion constant in Hz (set to 0 for simplicity)
        γX = -13.41679e6         # Spin-rotation interaction constant in Hz
        bFX = 170.26374e6        # Fermi contact term in Hz
        cX = 85.4028e6           # Dipolar hyperfine interaction constant in Hz
    end
    X_state_ham = Hamiltonian(basis=X_state_basis, operator=X_state_operator, parameters=X_state_parameters)
    Zeeman_x(state, state′) = (Zeeman(state, state′, -1) - Zeeman(state, state′, 1)) / sqrt(2)
    Zeeman_y(state, state′) = im * (Zeeman(state, state′, -1) + Zeeman(state, state′, 1)) / sqrt(2)
    Zeeman_z(state, state′) = Zeeman(state, state′, 0)

    X_state_ham = add_to_H(X_state_ham, :B_x, (gS * μB * 1e-6) * Zeeman_x)
    X_state_ham = add_to_H(X_state_ham, :B_y, (gS * μB * 1e-6) * Zeeman_y)
    X_state_ham = add_to_H(X_state_ham, :B_z, (gS * μB * 1e-6) * Zeeman_z)
    X_state_ham.parameters.B_x = 0.0
    X_state_ham.parameters.B_y = (1 / sqrt(2)) * B
    X_state_ham.parameters.B_z = (1 / sqrt(2)) * B

    # Evaluate and Solve the Hamiltonian
    # Compute the energy levels and eigenstates of the Hamiltonian.
    evaluate!(X_state_ham)
    QuantumStates.solve!(X_state_ham)
    all_grounds = X_state_ham.states


    # #################################################################################################
    # Excited state Hamiltonian
    # #################################################################################################
    QN_bounds = (
        S=1 / 2,            # Electron spin quantum number
        I=1 / 2,            # Nuclear spin quantum number
        Λ=(-1, 1),         # Projection of electronic orbital angular momentum (with possible values -1 and 1)
        J=1/2:3/2         # Total angular momentum quantum number range
    )

    # Enumerate states for the excited state basis
    A_state_basis = enumerate_states(HundsCaseA_LinearMolecule, QN_bounds)

    # Define the Hamiltonian operator for the excited state
    A_state_operator = :(
        T_A * DiagonalOperator +     # Diagonal constant (electronic zero point energy)
        Be_A * Rotation +            # Rotational term
        Aso_A * SpinOrbit +          # Spin-orbit interaction term
        q_A * ΛDoubling_q +          # Λ-doubling term q
        p_A * ΛDoubling_p2q +        # Λ-doubling term p2q
        a * Hyperfine_IL +           # Hyperfine interaction (IL term)
        d * Hyperfine_Dipolar_d      # Hyperfine interaction (dipolar term)
    )

    # Define parameters for the excited state Hamiltonian
    # Spectroscopic constants for CaOH, A state
    A_state_parameters = QuantumStates.@params begin
        T_A = 542.8102 * 10^12 # Diagonal constant (electrion zero point energy)
        Be_A = 40.9304 * 10^12 # Rotational constant
        Aso_A = 7.4276 * 10^9  # A spin-orbit constant
        p_A = 11.882 * 10^9
        q_A = 0
        a = -0.8e6
        d = -4.6e6
    end

    A_state_ham = Hamiltonian(basis=A_state_basis, operator=A_state_operator, parameters=A_state_parameters)
    evaluate!(A_state_ham)
    QuantumStates.solve!(A_state_ham)
    A_state_J12_pos_parity_states = A_state_ham.states[5:8]

    # Define quantum number bounds for the A state in Hund's case (b)
    QN_bounds = (
        S=1 / 2,
        I=1 / 2,
        Λ=(-1, 1),
        N=0:3 # ?
    )
    A_state_caseB_basis = enumerate_states(HundsCaseB_LinearMolecule, QN_bounds)

    ground_states = X_state_ham.states[5:16] # This is dropping all the terms with N = 0
    excited_states = convert_basis(A_state_J12_pos_parity_states, A_state_caseB_basis)

    states = [ground_states; excited_states]

    return ground_states, excited_states

end
