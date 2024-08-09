using Pkg
Pkg.activate(".")
using Revise, QuantumStates, OpticalBlochEquations, UnitsToValue, DifferentialEquations, Plots, PlotlyJS, Statistics, LinearAlgebra, StaticArrays, DataFrames, CSV, ProgressMeter
using StaticArrays, RectiGrids, StatsBase
RESULTS_PATH = "RESULTS"

#################################################################################################
# ## Helper Functions and Constants

# Function to get the max coefficient basis element
function get_max_coeff_basis_element(state)
    coeffs = abs.(state.coeffs) .^ 2
    max_idx = argmax(coeffs)
    return state.basis[max_idx]
end

# Function to extract quantum numbers from a state
function extract_quantum_numbers_from_state(state)
    basis_elem = get_max_coeff_basis_element(state)
    quantum_numbers = Dict(:N => basis_elem.N, :J => basis_elem.J, :F => basis_elem.F, :M => basis_elem.M)
    return quantum_numbers
end

# Function to find a state by quantum numbers
function find_state_by_quantum_numbers(states, N, J, F, M)
    for state in states
        qn = extract_quantum_numbers_from_state(state)
        if qn[:N] == N && qn[:J] == J && qn[:F] == F && qn[:M] == M
            return state
        end
    end
    error("State with quantum numbers N=$N, J=$J, F=$F, M=$M not found.")
end

# Function to get all state information
function all_state_info(state)
    basis_elem = get_max_coeff_basis_element(state)
    return basis_elem
end

_μB = (μB / h) * 1e-4
λ = 552e-9 # Wavelength of light in meters
Γ = 2π * 5.7e6
m = @with_unit 191 "u" # Mass of the molecule in atomic mass units
k = 2π / λ

#################################################################################################
# ## Setting up the states

# Define quantum number bounds for X and A states
function define_QN_bounds()
    QN_bounds_X = (
        S=1 / 2,
        I=1 / 2,
        Λ=0,
        N=0:3
    )

    QN_bounds_A = (
        S=1 / 2,
        I=1 / 2,
        Λ=(-1, 1),
        J=1/2:3/2
    )

    QN_bounds_A2 = (
        S=1 / 2,
        I=1 / 2,
        Λ=(-1, 1),
        N=0:3
    )

    return QN_bounds_X, QN_bounds_A, QN_bounds_A2
end

# Generate Hamiltonians for X and A states
function generate_hamiltonians(B)
    QN_bounds_X, QN_bounds_A, QN_bounds_A2 = define_QN_bounds()

    X_state_basis = enumerate_states(HundsCaseB_LinearMolecule, QN_bounds_X)
    X_state_operator = :(
        BX * Rotation +
        DX * RotationDistortion +
        γX * SpinRotation +
        bFX * Hyperfine_IS +
        cX * (Hyperfine_Dipolar / 3)
    )

    X_state_parameters = QuantumStates.@params begin
        BX = 7233.8271e6
        DX = 0.0
        γX = -13.41679e6
        bFX = 170.26374e6
        cX = 85.4028e6
    end

    X_state_ham = Hamiltonian(basis=X_state_basis, operator=X_state_operator, parameters=X_state_parameters)

    Zeeman_x(state, state′) = (Zeeman(state, state′, -1) - Zeeman(state, state′, 1)) / sqrt(2)
    Zeeman_y(state, state′) = im * (Zeeman(state, state′, -1) + Zeeman(state, state′, 1)) / sqrt(2)
    Zeeman_z(state, state′) = Zeeman(state, state′, 0)

    X_state_ham = add_to_H(X_state_ham, :B_x, (gS * _μB * 1e-6) * Zeeman_x)
    X_state_ham = add_to_H(X_state_ham, :B_y, (gS * _μB * 1e-6) * Zeeman_y)
    X_state_ham = add_to_H(X_state_ham, :B_z, (gS * _μB * 1e-6) * Zeeman_z)
    X_state_ham.parameters.B_x = 0.0
    X_state_ham.parameters.B_y = B / sqrt(2)
    X_state_ham.parameters.B_z = B / sqrt(2)

    evaluate!(X_state_ham)
    QuantumStates.solve!(X_state_ham)
    all_grounds = X_state_ham.states[5:16]
    ground_states = [state for state in all_grounds if extract_quantum_numbers_from_state(state)[:N] == 1]

    A_state_basis = enumerate_states(HundsCaseA_LinearMolecule, QN_bounds_A)

    A_state_operator = :(
        T_A * DiagonalOperator +
        Be_A * Rotation +
        Aso_A * SpinOrbit +
        q_A * ΛDoubling_q +
        p_A * ΛDoubling_p2q +
        a * Hyperfine_IL +
        d * Hyperfine_Dipolar_d
    )

    A_state_parameters = QuantumStates.@params begin
        T_A = 542.8102 * 10^12
        Be_A = 40.9304 * 10^12
        Aso_A = 7.4276 * 10^9
        p_A = 11.882 * 10^9
        q_A = 0
        a = -0.8e6
        d = -4.6e6
    end

    A_state_ham = Hamiltonian(basis=A_state_basis, operator=A_state_operator, parameters=A_state_parameters)
    evaluate!(A_state_ham)
    QuantumStates.solve!(A_state_ham)

    excited_states = A_state_ham.states[1:4]
    A_state_basis = enumerate_states(HundsCaseB_LinearMolecule, QN_bounds_A2)
    excited_states = convert_basis(excited_states, A_state_basis)

    return ground_states, excited_states, X_state_ham, A_state_ham, A_state_basis
end

# Generate Transition Dipole Moments (TDM)
function generate_tdms(X_state_ham, ground_states, excited_states, A_state_basis)
    d = zeros(ComplexF64, 16, 16, 3)
    d_ge = zeros(ComplexF64, 12, 4, 3)
    basis_tdms = get_tdms_two_bases(X_state_ham.basis, A_state_basis, TDM)
    tdms_between_states!(d_ge, basis_tdms, ground_states, excited_states)
    d[1:12, 13:16, :] .= d_ge
    return d
end


#################################################################################################
# ## Setting up the lasers

# Define laser parameters and generate lasers
function generate_lasers(config, pol_J12, ω_J12s, s_J12_sidebands)
    ϵ(ϵ_val) = t -> ϵ_val
    use_single_state = false
    polarizations = [pol_J12 for _ in ω_J12s]
    if use_single_state

        # Define the wave vectors and polarizations for the lasers in the XY‖ configuration
        k̂ = +x̂
        ϵ1 = ϵ(rotate_pol(pol_J12, ẑ))
        laser1 = Field(k̂, ϵ1, ω_J12, s_J12)
        k̂ = -x̂
        ϵ2 = ϵ(rotate_pol(pol_J12, -ẑ))
        laser2 = Field(k̂, ϵ2, ω_J12, s_J12)
        k̂ = +ŷ
        ϵ3 = ϵ(rotate_pol(pol_J12, ẑ))
        laser3 = Field(k̂, ϵ3, ω_J12, s_J12)
        k̂ = -ŷ
        ϵ4 = ϵ(rotate_pol(pol_J12, -ẑ))
        laser4 = Field(k̂, ϵ4, ω_J12, s_J12)

        lasers_XY_parallel = [laser1, laser2, laser3, laser4]
        # Define the wave vectors and polarizations for the lasers in the XY⊥ configuration
        k̂ = +x̂
        ϵ1 = ϵ(rotate_pol(pol_J12, ŷ))
        laser1 = Field(k̂, ϵ1, ω_J12, s_J12)
        k̂ = -x̂
        ϵ2 = ϵ(rotate_pol(pol_J12, -ŷ))
        laser2 = Field(k̂, ϵ2, ω_J12, s_J12)
        k̂ = +ŷ
        ϵ3 = ϵ(rotate_pol(pol_J12, ẑ))
        laser3 = Field(k̂, ϵ3, ω_J12, s_J12)
        k̂ = -ŷ
        ϵ4 = ϵ(rotate_pol(pol_J12, -ẑ))
        laser4 = Field(k̂, ϵ4, ω_J12, s_J12)

        lasers_XY_perpendicular = [laser1, laser2, laser3, laser4]
    else

        lasers_plus_x = [Field(+x̂, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, ẑ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        lasers_minus_x = [Field(-x̂, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, -ẑ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        lasers_plus_y = [Field(+ŷ, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, ẑ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        lasers_minus_y = [Field(-ŷ, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, -ẑ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        lasers_XY_parallel = [lasers_plus_x; lasers_minus_x; lasers_plus_y; lasers_minus_y]

        laser_plus_x = [Field(+x̂, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, ŷ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        laser_minus_x = [Field(-x̂, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, -ŷ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        laser_plus_y = [Field(+ŷ, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, ẑ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        laser_minus_y = [Field(-ŷ, ϵ, ω, s) for (ϵ, j, ω, s) in zip([ϵ(rotate_pol(pol_J12, -ẑ)) for j in 1:length(ω_J12s)], polarizations, ω_J12s, s_J12_sidebands)]
        lasers_XY_perpendicular = [laser_plus_x; laser_minus_x; laser_plus_y; laser_minus_y]
    end
    return lasers_XY_parallel, lasers_XY_perpendicular
end

#################################################################################################
# Define and solve Optical Bloch Equations (OBE)
function solve_obe(states, lasers, d)
    ρ0 = zeros(ComplexF64, length(states), length(states))
    particle = OpticalBlochEquations.Particle()

    ρ0[1, 1] = 1.0

    function prob_func!(prob, scan_values_grid, i)
        p = prob.p
        p.v .= (0, scan_values_grid[i].v, 0)
        p.v .= round_vel(p.v, p.freq_res)
        p.r0 .= scan_values_grid[i].r
        return prob
    end
    function output_func(p, sol)
        f = p.force_last_period
        return (f[1], f[2], f[3])
    end
    freq_res = 1e-1
    p = obe(ρ0, particle, states, lasers, d, true, true; λ=λ, Γ=Γ, freq_res=freq_res)

    t_end = 5p.period + 1
    tspan = (0.0, t_end)
    @time prob = ODEProblem(ρ!, p.ρ0_vec, tspan, p, reltol=1e-3, save_on=false)

    di = 2
    rs = vcat([(n1 / (di + 1), n2 / (di + 1), n3 / (di + 1)) .* 2π for n1 ∈ 0:di, n2 ∈ 0:di, n3 ∈ 0:di]...)
    vs = 0.0:0.5:5

    scan_values = (r=rs, v=vs)
    scan_values_grid = RectiGrids.grid(scan_values)
    forces, populations = force_scan_v2(prob, scan_values_grid, prob_func!, output_func)
    return forces, populations
end

function get_omegas(ground_states, excited_states, detunning, Γ)
    # source states
    g_interest_states = [state for state in ground_states if (all_state_info(state).v_1 == 0 && all_state_info(state).v_2 == 0 && all_state_info(state).v_3 == 0)][1:4]
    input_energies = [energy(st) for st in g_interest_states]
    # Target states
    a_qns = [extract_quantum_numbers_from_state(st) for st in excited_states][1]
    energies = [energy(st) for st in excited_states]
    # sorted
    target_states = sort(excited_states, by=energy)
    target_energy = [energy(st) for st in excited_states][1]
    A_energy = target_energy
    δJ12 = detunning * Γ
    ω_J12s = [2π * (A_energy - inp_e) + δJ12 for inp_e in input_energies]
    return ω_J12s
end

function save_as_df(forces, populations, saving_path)
    averaged_forces = []
    for (i, v) ∈ enumerate(vs)
        idxs = [j for (j, x) ∈ enumerate(scan_values_grid) if x.v == v]
        _forces = [f[1] for f in forces[idxs] if abs(f[1]) <= 1e3]
        push!(averaged_forces, mean([f[1] for f in _forces]))
    end
    df = DataFrame(v=vs, force=averaged_forces, population=populations)
    CSV.write(saving_path, df)
end

# Main function to run the simulation
function run_simulation(detunings, magnetic_fields, CONFIGURATION)
    Γ = 2π * 5.7e6
    for B in magnetic_fields
        for detuning in detunings
            @show B, detuning
            ground_states, excited_states, X_state_ham, A_state_ham, A_state_basis = generate_hamiltonians(10.0 * 1e-4)
            d = generate_tdms(X_state_ham, ground_states, excited_states, A_state_basis)
            ω_J12s = get_omegas(ground_states, excited_states, detuning, Γ)
            s_func(s) = (r, t) -> s
            I_sat = π * h * c * Γ / (3λ^3) / 10
            s_J12 = s_func(100.0)
            ratios = [1, 1 / 2, 1 / 2, 0]
            s_J12_sidebands = [s_func(ratio * 100.0) for ratio in ratios]
            lasers_XY_parallel, lasers_XY_perpendicular = generate_lasers(CONFIGURATION, σ⁺, ω_J12s, s_J12_sidebands)
            if CONFIGURATION == "XY_parallel"
                lasers = lasers_XY_parallel
            else
                lasers = lasers_XY_perpendicular
            end
            states = [ground_states; excited_states]
            forces, populations = solve_obe(states, lasers, d)
            save_as_df(forces, populations, "data/$(CONFIGURATION)_B_$(B)_detuning_$(detuning).csv")
        end
    end
end


detunings = [1.0, 3.0, 5.0, 10.0] * 1e6
magnetic_fields = [0.0, 1.0, 5.0, 25.0] .* 1e-4
magnetic_fields = magnetic_fields[2:2]
run_simulation(detunings, magnetic_fields, "XY_parallel")
