using Statistics, LinearAlgebra

number_of_particles = 100000  # Number of particles to input

input_particle = "e-"       # e- = electrons, proton = protons, gamma = photons

energy_kev_min = 10         # Minimum beam energy, keV
energy_kev_max = 10_000     # Maximum beam energy, keV
energy_nbeams = 5          # Number of log-spaced beams to place between minimum and maximum energy

pitch_angle_deg_min = 0     # Minimum beam pitch angle, deg
pitch_angle_deg_max = 90    # Maximum beam pitch angle, deg
pitch_angle_nbeams = 6     # Number of linear-spaced beams to place between minimum and maximum pitch angle

# Create energy and pitch angle lists
energies_to_simulate = 10.0 .^ LinRange(log10(energy_kev_min), log10(energy_kev_max), energy_nbeams)
pitch_angles_to_simulate = LinRange(pitch_angle_deg_min, pitch_angle_deg_max, pitch_angle_nbeams)

# Round energies and pitch angles to nearest integer for reduced filename verbosity. If you have sub-integer resolution,
# this will eat that, so if you're doing < 1º pitch angle resolution or < 1 keV energy resolution, remove this block. Note
# that removal of this block may result in very long output filenames.
energies_to_simulate = round.(energies_to_simulate)
pitch_angles_to_simulate = round.(pitch_angles_to_simulate)

# Create shell script
file = open("$(@__DIR__)/all_beams.txt", "w")
[println(file, "./G4EPP $(number_of_particles) $(input_particle) $(E) $(α)") for E in energies_to_simulate, α in pitch_angles_to_simulate]
close(file)