"""
Function to calculate detuning values from experimental voltage and AOM frequency.
:param voltage: Voltage value.
:param aom_freq: AOM frequency.
:return: Δ1, Δ2 in MHz.
"""
function get_Δ_from_exp(voltage, aom_freq)
    Δ1 = 57 - 7.4 * (5.5 - voltage)
    Δ2 = Δ1 + 51.24 - aom_freq
    return Δ1, Δ2
end

"""
Function to flip the polarization vector components.
:param ϵ: Polarization vector.
:return: Flipped polarization vector.
"""
function flip(ϵ)
    return SVector{3,ComplexF64}(ϵ[3], ϵ[2], ϵ[1])
end

"""
Function to calculate Gaussian intensity along specified axes.
:param r: Position vector.
:param axes: Axes to calculate intensity along.
:param centers: Centers of the Gaussian beams.
:return: Gaussian intensity value.
"""
function gaussian_intensity_along_axes(r, axes, centers)
    """1/e^2 width = 5mm Gaussian beam """
    d2 = (r[axes[1]] - centers[1])^2 + (r[axes[2]] - centers[2])^2
    return exp(-2 * d2 / (5e-3 / (1 / k))^2)
end

"""
Function to define lasers with specific parameters and configurations.
:param states: Energy states.
:param RF_frequency: RF frequency for modulation.
:param s1, s2, s3, s4: Laser intensities.
:param Δ1, Δ2, Δ3, Δ4: Detuning values.
:param pol_J12: Initial polarization for J=1/2 state.
:param pol_imbalance: Polarization imbalance factor.
:param s_imbalance: Laser intensity imbalance factor.
:param retro_loss: Retroreflection loss factor.
:param off_center: Off-center displacement values.
:param pointing_error: Pointing error values.
:return: Two lists of laser configurations: XY parallel and XY perpendicular.
"""
function define_lasers(
    states,
    RF_frequency,
    s1,
    s2,
    s3,
    s4,
    Δ1,
    Δ2,
    Δ3,
    Δ4,
    pol_J12, # New parameter for defining the initial polarization
    pol_imbalance,
    s_imbalance,
    retro_loss,
    off_center,
    pointing_error
)

    # Calculate angular frequencies for each laser
    ω1 = 2π * (energy(states[end]) - energy(states[1])) + Δ1
    ω2 = 2π * (energy(states[end]) - energy(states[5])) + Δ2
    ω3 = 2π * (energy(states[end]) - energy(states[5])) + Δ3
    ω4 = 2π * (energy(states[end]) - energy(states[5])) + Δ4

    # Random centers for Gaussian beams
    x_center_y = rand() * off_center[1] * k
    x_center_z = rand() * off_center[2] * k
    y_center_x = rand() * off_center[3] * k
    y_center_z = rand() * off_center[4] * k
    z_center_x = rand() * off_center[5] * k
    z_center_y = rand() * off_center[6] * k

    # Function to modulate polarization based on time and RF frequency
    ϵ_(ϵ1, ϵ2) = t -> iseven(t ÷ (π / RF_frequency)) ? ϵ1 : ϵ2

    # Gaussian intensity function
    s_gaussian(s, axes, centers) = (r, t) -> s * gaussian_intensity_along_axes(r, axes, centers)

    # Adjust polarizations for imbalance
    rand1 = rand()
    pol_J12_x = pol_J12 .* sqrt(1 - rand1 * pol_imbalance) + flip(pol_J12) .* sqrt(rand1 * pol_imbalance)
    rand2 = rand()
    pol_J12_y = pol_J12 .* sqrt(1 - rand2 * pol_imbalance) + flip(pol_J12) .* sqrt(rand2 * pol_imbalance)

    # Random intensity imbalances
    sx_rand = 1 / 2 - rand()
    sy_rand = 1 / 2 - rand()
    sz_rand = 1 / 2 - rand()

    # Random phase shifts
    ϕs = [exp(im * 2π * rand()), exp(im * 2π * rand()), exp(im * 2π * rand()), exp(im * 2π * rand()), exp(im * 2π * rand()), exp(im * 2π * rand())]

    # Wave vectors with pointing errors
    kx = x̂ + [0, pointing_error[1], pointing_error[2]]
    kx = kx ./ sqrt(kx[1]^2 + kx[2]^2 + kx[3]^2)
    ky = ŷ + [pointing_error[3], 0.0, pointing_error[4]]
    ky = ky ./ sqrt(ky[1]^2 + ky[2]^2 + ky[3]^2)
    kz = ẑ + [pointing_error[5], pointing_error[6], 0.0]
    kz = kz / sqrt(kz[1]^2 + kz[2]^2 + kz[3]^2)

    # Intensity adjustments for each laser
    s1x = s1 * (1 + s_imbalance[1] * sx_rand)
    s1y = s1 * (1 + s_imbalance[2] * sy_rand)
    s1z = s1 * (1 + s_imbalance[3] * sz_rand)

    # Polarization functions for XY‖ configuration
    ϵ_func_x = ϵ_(rotate_pol(pol_J12_x, ẑ), rotate_pol(pol_J12_x, -ẑ))
    ϵ_func_y = ϵ_(rotate_pol(pol_J12_y, ẑ), rotate_pol(pol_J12_y, -ẑ))

    # Define lasers for XY‖ configuration
    k̂1 = +kx
    laser1 = Field(k̂1, ϵ_func_x, ω1, s_gaussian(s1x, (2, 3), (x_center_y, x_center_z)))
    k̂2 = -kx
    laser2 = Field(k̂2, ϵ_func_x, ω1, s_gaussian(s1x * (1 - retro_loss), (2, 3), (x_center_y, x_center_z)))
    k̂3 = +ky
    laser3 = Field(k̂3, ϵ_func_y, ω1, s_gaussian(s1y, (1, 3), (y_center_x, y_center_z)))
    k̂4 = -ky
    laser4 = Field(k̂4, ϵ_func_y, ω1, s_gaussian(s1y * (1 - retro_loss), (1, 3), (y_center_x, y_center_z)))

    lasers_XY_parallel = [laser1, laser2, laser3, laser4]

    # Polarization functions for XY⊥ configuration
    ϵ_func_x_perp = ϵ_(rotate_pol(pol_J12_x, ŷ), rotate_pol(pol_J12_x, -ŷ))
    ϵ_func_y_perp = ϵ_(rotate_pol(pol_J12_y, -ẑ), rotate_pol(pol_J12_y, ẑ))

    # Define lasers for XY⊥ configuration
    k̂1 = +kx
    laser1 = Field(k̂1, ϵ_func_x_perp, ω1, s_gaussian(s1x, (2, 3), (x_center_y, x_center_z)))
    k̂2 = -kx
    laser2 = Field(k̂2, ϵ_func_x_perp, ω1, s_gaussian(s1x * (1 - retro_loss), (2, 3), (x_center_y, x_center_z)))
    k̂3 = +ky
    laser3 = Field(k̂3, ϵ_func_y_perp, ω1, s_gaussian(s1y, (1, 3), (y_center_x, y_center_z)))
    k̂4 = -ky
    laser4 = Field(k̂4, ϵ_func_y_perp, ω1, s_gaussian(s1y * (1 - retro_loss), (1, 3), (y_center_x, y_center_z)))

    lasers_XY_perpendicular = [laser1, laser2, laser3, laser4]

    return lasers_XY_parallel, lasers_XY_perpendicular
end
