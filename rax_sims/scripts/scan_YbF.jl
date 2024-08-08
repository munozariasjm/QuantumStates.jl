# %% [markdown]
# # YbF
#
# This script aims to study the effects of detuning and magnetic field on the force profile of YbF molecules, based on the results from [An ultracold molecular beam for testing fundamental physics](https://iopscience.iop.org/article/10.1088/2058-9565/ac107e).

# %% [markdown]
# ## Loading packages

# %%
using Pkg
Pkg.activate("/Users/jose/Documents/Works/MIT/RaX/Simu/Molecule-Sims/")
using Revise, QuantumStates, OpticalBlochEquations, UnitsToValue, DifferentialEquations, Plots, PlotlyJS, Statistics, LinearAlgebra, StaticArrays, DataFrames

# %% [markdown]
# ## Helper Functions and Constants

# %%
function get_max_coeff_basis_element(state)
    coeffs = abs.(state.coeffs) .^ 2
    max_idx = argmax(coeffs)
    return state.basis[max_idx]
end

function extract_quantum_numbers_from_state(state)
    basis_elem = get_max_coeff_basis_element(state)
    quantum_numbers = Dict(:N => basis_elem.N, :J => basis_elem.J, :F => basis_elem.F, :M => basis_elem.M)
    return quantum_numbers
end

function find_state_by_quantum_numbers(states, N, J, F, M)
    for state in states
        qn = extract_quantum_numbers_from_state(state)
        if qn[:N] == N && qn[:J] == J && qn[:F] == F && qn[:M] == M
            return state
        end
    end
    error("State with quantum numbers N=$N, J=$J, F=$F, M=$M not found.")
end

function all_state_info(state)
    basis_elem = get_max_coeff_basis_element(state)
    return basis_elem
end

_μB = (μB / h) * 1e-4;

# %% [markdown]
# ## Problem Setup

# %%
# Quantum Number Bounds
QN_bounds = (
    S=1 / 2,
    I=1 / 2,
    Λ=0,
    N=0:3
)

X_state_basis = enumerate_states(HundsCaseB_LinearMolecule, QN_bounds)
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

B = 10.0 * 1e-4 # 10 Gau

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

QN_bounds = (
    S=1 / 2,
    I=1 / 2,
    Λ=(-1, 1),
    J=1/2:3/2
)
A_state_basis = enumerate_states(HundsCaseA_LinearMolecule, QN_bounds)

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
QN_bounds = (
    S=1 / 2,
    I=1 / 2,
    Λ=(-1, 1),
    N=0:3
)
A_state_basis = enumerate_states(HundsCaseB_LinearMolecule, QN_bounds)
excited_states = convert_basis(excited_states, A_state_basis)
[extract_quantum_numbers_from_state(state) for state in excited_states]

d = zeros(ComplexF64, 16, 16, 3)
d_ge = zeros(ComplexF64, 12, 4, 3)
basis_tdms = get_tdms_two_bases(X_state_ham.basis, A_state_basis, TDM)
tdms_between_states!(d_ge, basis_tdms, ground_states, excited_states)
d[1:12, 13:16, :] .= d_ge

states = [ground_states; excited_states]
n_excited = length(excited_states)

λ = 552e-9
Γ = 2π * 5.7e6
m = @with_unit 191 "u"
k = 2π / λ
δJ12 = +6 * Γ

g_interest_states = [state for state in all_grounds if (all_state_info(state).v_1 == 0 && all_state_info(state).v_2 == 0 && all_state_info(state).v_3 == 0)][1:4]
input_energies = [energy(st) for st in g_interest_states]

a_qns = [extract_quantum_numbers_from_state(st) for st in excited_states][1]
energies = [energy(st) for st in excited_states]
target_states = sort(excited_states, by=energy)
target_energy = [energy(st) for st in excited_states][1]

delta_state_energies = target_energy .- input_energies
ω_J12s = [2π * (target_energy - inp_e) + δJ12 for inp_e in input_energies]
ω_J12 = ω_J12s[1]
s_func(s) = (r, t) -> s
I_sat = π * h * c * Γ / (3λ^3) / 10
s_J12 = s_func(100.0)
ratios = [1, 1 / 2, 1 / 2, 0]
s_J12_sidebands = [s_func(ratio * 100.0) for ratio in ratios]

pol_J12 = σ⁺
polarizations = [pol_J12 for _ in 1:length(ω_J12s)]
ϵ(ϵ_val) = t -> ϵ_val;

k̂ = +x̂;
ϵ1 = ϵ(rotate_pol(pol_J12, ẑ));
laser1 = Field(k̂, ϵ1, ω_J12, s_J12);
k̂ = -x̂;
ϵ2 = ϵ(rotate_pol(pol_J12, -ẑ));
laser2 = Field(k̂, ϵ2, ω_J12, s_J12);
k̂ = +ŷ;
ϵ3 = ϵ(rotate_pol(pol_J12, ẑ));
laser3 = Field(k̂, ϵ3, ω_J12, s_J12);
k̂ = -ŷ;
ϵ4 = ϵ(rotate_pol(pol_J12, -ẑ));
laser4 = Field(k̂, ϵ4, ω_J12, s_J12);

lasers_XY_parallel = [laser1, laser2, laser3, laser4]

lasers = lasers_XY_parallel

particle = OpticalBlochEquations.Particle()
ρ0 = zeros(ComplexF64, length(states), length(states))
names = [extract_quantum_numbers_from_state(st) for st in states]
string_names = ["|$(qn[:J]),  $(qn[:F]), $(qn[:M])>" for qn in names]
ρ0[1, 1] = 1.0

freq_res = 1e-2
p = obe(ρ0, particle, states, lasers, d, true, true; λ=λ, Γ=Γ, freq_res=freq_res)

bounds_FWHMs = (1.5, 3.1)
mean_fw = mean(bounds_FWHMs)
function fwhm_to_rms(x)
    return x / (2 * sqrt(2 * log(2)))
end
rms_speed = fwhm_to_rms(mean_fw)

p.r0 = (0.0, 0.0, 0.0)
p.v = (0.0, 0.0, rms_speed) ./ (Γ / k)
p.v = round_vel(p.v, p.freq_res)

t_end = 5p.period + 1
tspan = (0.0, t_end)
prob = ODEProblem(ρ!, p.ρ0_vec, tspan, p)
times = range(0, t_end, 100_00)

cb = PeriodicCallback(reset_force!, p.period)
@time sol = DifferentialEquations.solve(prob, DP5(), callback=cb, reltol=1e-3, saveat=times)

function plot_force_profile(detuning, B)
    p.δ = detuning
    X_state_ham.parameters.B_y = B / sqrt(2)
    X_state_ham.parameters.B_z = B / sqrt(2)

    prob = ODEProblem(ρ!, p.ρ0_vec, tspan, p)
    sol = DifferentialEquations.solve(prob, DP5(), callback=cb, reltol=1e-3, saveat=times)

    plot_us = sol.u
    plot_ts = sol.t
    n_states = size(p.ρ_soa, 1)

    Plots.plot(
        ylim=(-0.0, 1.6),
    )
    for i in 1:n_states
        state_idx = n_states * (i - 1) + i
        if maximum([real(u[state_idx]) for u in plot_us]) > 0.1
            Plots.plot!(plot_ts, [real(u[state_idx]) for u in plot_us], label=string_names[i])
        else
            Plots.plot!(plot_ts, [real(u[state_idx]) for u in plot_us], label=false)
        end
    end
    Plots.plot!(
        xlabel="Time (s)",
        ylabel="Population",
        title="Force Profile (Detuning: $detuning, B-field: $B)"
    )
end

detunings = [-6 * Γ, -3 * Γ, 0, 3 * Γ, 6 * Γ]
magnetic_fields = [5.0 * 1e-4, 10.0 * 1e-4, 15.0 * 1e-4]

for detuning in detunings
    for B in magnetic_fields
        plot_force_profile(detuning, B)
    end
end
