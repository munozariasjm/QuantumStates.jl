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