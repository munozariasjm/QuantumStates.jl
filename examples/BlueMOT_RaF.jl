using Pkg
Pkg.activate("/Users/jose/Documents/Works/MIT/RaX/Simu/Molecule-Sims/")

# Import the primary packages required for the notebook
using
    Revise,
    QuantumStates,            # for calculating molecular structure
    OpticalBlochEquations,    # for solving optical Bloch equations
    UnitsToValue              # for numerical values
;

### Create Hamiltonian for the $X^2\Sigma^+(000, N=1)$ ground state

# We first create a Hamiltonian for the $N=1$ state of the ground electronic $X^2\Sigma^+(000)$ state. Creating a Hamiltonian requiries defining three distinct objects:
#     1) A basis, accomplished using the function `enumerate_states`, which takes a basis state type as its first argument and appropriate bounds for the quantum numbers. Note that we let the basis include states of $N \in [0,3]$, not just $N=1$, since we could have matrix elements connecting the $N=1$ states to other rotational states. Such matrix elements are only included in the Hamiltonian provided that the basis includes all associated basis states.
#     2) An operator for the Hamiltonian based on matrix elements that have been defined for the given basis state type. The syntax is to write the operator as a sum of terms of the form [parameter symbol] * [matrix operator], e.g., `BX * Rotation`.
#     3) A set of parameters using the macro `QuantumStates.@params` (note that other packages export a `@params` macro, hence the inclusion of the package name; normally we could just write `@params`). These parameters must correspond one-to-one with the parameter symbols used in the definition of the operator.

#     Additional terms can later be added to the Hamiltonian object using the `add_to_H` function, see below for an example.

#     Finally, the Hamiltonian is evaluated using `evaluate!` and then solved using `QuantumStates.solve!` (other packages export a `solve!` function, similarly to `@params` and the package name is therefore included in the function call).

QN_bounds = (
    S=1 / 2,
    I=1 / 2,
    Λ=0,
    N=0:3
)
X_state_basis = enumerate_states(HundsCaseB_LinearMolecule, QN_bounds)

X_state_operator = :(
    BX * Rotation +                     #
    DX * RotationDistortion +
    γX * SpinRotation #+                 # Spin-rotation interaction
    # bFX * Hyperfine_IS +
    # cX * (Hyperfine_Dipolar/3)
)

# Now we have to plug in the constants for RaF:
# Taken from: https://www.nature.com/articles/s41567-023-02296-w

X_state_parameters = QuantumStates.@params begin
    BX = 0.191985
    DX = 1.405 * (1e-7)
    γX = 0.00585
    bFX = 0
    cX = 0
end

X_state_ham = Hamiltonian(basis=X_state_basis, operator=X_state_operator, parameters=X_state_parameters)

# Add Zeeman terms
Zeeman_z(state, state′) = Zeeman(state, state′, 0)
X_state_ham = add_to_H(X_state_ham, :B_z, (gS * μB / h) * Zeeman_z)
X_state_ham.parameters.B_z = 0.0

evaluate!(X_state_ham)
QuantumStates.solve!(X_state_ham)
;