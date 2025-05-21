using Statistics, LinearAlgebra

number_of_particles = Int(1e5)  # Number of particles to input

input_particle = "e-"       # e- = electrons, proton = protons, gamma = photons

energy_kev_min = 10         # Minimum beam energy, keV
energy_kev_max = 10_000     # Maximum beam energy, keV
energy_nbeams = 5           # Number of log-spaced beams to place between minimum and maximum energy

pitch_angle_deg_min = 0     # Minimum beam pitch angle, deg
pitch_angle_deg_max = 90    # Maximum beam pitch angle, deg
pitch_angle_nbeams = 9      # Number of linear-spaced beams to place between minimum and maximum pitch angle

# Create energy and pitch angle lists
energies_to_simulate = 10.0 .^ LinRange(log10(energy_kev_min), log10(energy_kev_max), energy_nbeams)
pitch_angles_to_simulate = LinRange(pitch_angle_deg_min, pitch_angle_deg_max, pitch_angle_nbeams)

# Round energies and pitch angles to nearest integer for reduced filename verbosity. If you have sub-integer resolution,
# this will eat that, so if you're doing < 1º pitch angle resolution or < 1 keV energy resolution, remove this block. Note
# that removal of this block may result in very long output filenames.
energies_to_simulate = round.(energies_to_simulate)
pitch_angles_to_simulate = round.(pitch_angles_to_simulate)





# Create shell script
E = energies_to_simulate[1]
α = pitch_angles_to_simulate[1]

job_name = "$(E)keV_$(α)deg"

file = open("$(@__DIR__)/$(job_name).sh", "w")
println(file,
"""
#!/bin/bash

#SBATCH --job-name G4EPP_$(job_name)
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 05:00:00
#SBATCH --output /projects/jucl6426/G4EPP/results/batch/log_$(job_name).out
#SBATCH --qos=blanca-lair
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Run simulation
cd /projects/jucl6426/G4EPP/build/
./G4EPP $(number_of_particles) e- $(E) $(α)

# Copy results to safe folder
cp -r /projects/jucl6426/G4EPP/build/results/input_450.0km_record_450.0km/backscatter_electron_input_$(E)keV_$(α)deg_$(number_of_particles)particles /projects/jucl6426/G4EPP/results/batch
cp -r /projects/jucl6426/G4EPP/build/results/input_450.0km_record_450.0km/energy_deposition_electron_input_$(E)keV_$(α)deg_$(number_of_particles)particles /projects/jucl6426/G4EPP/results/batch
"""
)
close(file)