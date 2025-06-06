using Statistics, LinearAlgebra
using Glob

number_of_particles = 100_000  # Number of particles to input

input_particle = "e-"          # "e-" = electrons, "proton" = protons, "gamma" = photons

energy_kev_min = 10            # Minimum beam energy, keV
energy_kev_max = 10_000        # Maximum beam energy, keV
energy_nbeams = 10             # Number of log-spaced beams to place between minimum and maximum energy

pitch_angle_deg_min = 0        # Minimum beam pitch angle, deg
pitch_angle_deg_max = 90       # Maximum beam pitch angle, deg
pitch_angle_nbeams = 10        # Number of linear-spaced beams to place between minimum and maximum pitch angle

# Create energy and pitch angle lists
energies_to_simulate = 10.0 .^ LinRange(log10(energy_kev_min), log10(energy_kev_max), energy_nbeams)
pitch_angles_to_simulate = LinRange(pitch_angle_deg_min, pitch_angle_deg_max, pitch_angle_nbeams)

# TODO delete this. this is just for ELFIN purposes
pitch_angles_to_simulate = float.([0:5:60..., 61:1:69..., 70:5:90...])
energies_to_simulate = [63.245540618896484, 97.97958374023438, 138.5640869140625, 183.30308532714844, 238.11758422851562, 305.20489501953125, 385.16229248046875, 520.48046875, 752.9939575195312, 1081.665283203125, 1529.7060546875, 2121.3203125, 2893.960205078125, 3728.6064453125, 4906.12060546875, 6500.0]

# Round energies and pitch angles to nearest integer for reduced filename verbosity. If you have sub-integer resolution,
# this will eat that, so if you're doing < 1º pitch angle resolution or < 1 keV energy resolution, remove this block. Note
# that removal of this block may result in very long output filenames.
energies_to_simulate = round.(energies_to_simulate)
pitch_angles_to_simulate = round.(pitch_angles_to_simulate)

# Don't simulate if we already have data for a given beam
include("/Users/luna/Research/G4EPP_2.0/Frontend_Functions.jl")
energies_to_remove, pitch_angles_to_remove = _get_beam_locations(data_type = "raw")
skipped = 0

# Create shell scripts
rm.(glob("*deg.sh", @__DIR__))
for E in energies_to_simulate
  for α in pitch_angles_to_simulate
    if (E, α) ∈ collect(zip(energies_to_remove, pitch_angles_to_remove))
      global skipped += 1
      continue
    end

    job_name = "$(E)keV_$(α)deg"
    qos = "preemptable"
    time_limit = "1-00:00:00"

    if (α > 65) && (E < 400)
      qos = "blanca-lair"
      time_limit = "7-00:00:00"
    end

    file = open("$(@__DIR__)/$(job_name).sh", "w")
    println(file,
    """
    #!/bin/bash

    #SBATCH --job-name G4EPP_$(job_name)
    #SBATCH --nodes 1
    #SBATCH --ntasks-per-node 40
    #SBATCH --time $(time_limit)
    #SBATCH --output /projects/jucl6426/G4EPP/results/log_$(job_name).out
    #SBATCH --qos=$(qos)
    #SBATCH --exclude=bhpc-c5-u7-22,bhpc-c5-u7-23
    #SBATCH --no-requeue
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=jucl6426@colorado.edu

    # Terminate on any non-zero exit status
    set -e

    # Run simulation
    cd /projects/jucl6426/G4EPP/build/
    ./G4EPP $(number_of_particles) e- $(E) $(α)

    # Copy results to safe folder
    cp /projects/jucl6426/G4EPP/build/results/input_449.5km_record_450.5km/backscatter_electron_input_$(E)keV_$(α)deg_$(number_of_particles)particles.csv /projects/jucl6426/G4EPP/results
    cp /projects/jucl6426/G4EPP/build/results/input_449.5km_record_450.5km/energy_deposition_electron_input_$(E)keV_$(α)deg_$(number_of_particles)particles.csv /projects/jucl6426/G4EPP/results
    """
    )
    close(file)
  end
end

println("Skipped $(skipped)/$(length(energies_to_simulate)*length(pitch_angles_to_simulate)) files.")