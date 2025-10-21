using Statistics, LinearAlgebra
using Glob
using Printf
include("/Users/luna/Research/G4EPP_2.0/G4EPP.jl")

number_of_particles = 100_000  # Number of particles to input

particle = "e-"          # "e-" = electrons, "proton" = protons, "gamma" = photons

# Create energy and pitch angle lists
energy_kev_min = 1            # Minimum beam energy, keV
energy_kev_max = 10_000       # Maximum beam energy, keV
energy_nbeams = 30             # Number of log-spaced beams to place between minimum and maximum energy
energies_to_simulate = 10.0 .^ LinRange(log10(energy_kev_min), log10(energy_kev_max), energy_nbeams)

pitch_angle_deg_min = 0        # Minimum beam pitch angle, deg
pitch_angle_deg_max = 90       # Maximum beam pitch angle, deg
pitch_angle_nbeams = 19        # Number of linear-spaced beams to place between minimum and maximum pitch angle
pitch_angles_to_simulate = LinRange(pitch_angle_deg_min, pitch_angle_deg_max, pitch_angle_nbeams)


# Round energies and pitch angles to nearest integer for reduced filename verbosity. If you have sub-integer resolution,
# this will eat that, so if you're doing < 1º pitch angle resolution or < 1 keV energy resolution, remove this block. Note
# that removal of this block may result in very long output filenames.
energies_to_simulate = round.(energies_to_simulate, digits = 1)
pitch_angles_to_simulate = round.(pitch_angles_to_simulate, digits = 1)

# Find beams we already data for
prev_particle_type, prev_energy, prev_pitch_angle = get_beams("/Users/luna/Research/geant4/G4EPP/_Full_Grid")

# Create shell scripts
rm.(glob("*deg.sh", @__DIR__))
written = 0
skipped = 0

for E in energies_to_simulate
  for α in pitch_angles_to_simulate
    input_particle_longname = particle == "e-" ? "electron" : particle
    energy_string = @sprintf "%.1f" E
    job_name = "$(input_particle_longname)_$(energy_string)keV_$(α)deg"
    #qos = "preemptable"
    #time_limit = "1-00:00:00"
    qos = "blanca-lair"
    time_limit = "7-00:00:00"
    
    # Don't simulate if we already have data for a given beam
    if sum( (input_particle_longname .== prev_particle_type) .&& (E .== prev_energy) .&& (α .== prev_pitch_angle) ) == 1
      continue
    end

    if E > 500
      continue
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
    #SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-20,bhpc-c5-u7-22,bhpc-c5-u7-23
    #SBATCH --requeue
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=jucl6426@colorado.edu

    # Terminate on any non-zero exit status
    set -e

    # Run simulation
    cd /projects/jucl6426/G4EPP/build/
    ./G4EPP $(number_of_particles) $(particle) $(energy_string) $(α)

    # Copy results to safe folder
    cp /projects/jucl6426/G4EPP/build/results/mlat_65.77deg_input_449.5km_record_450.5km/backscatter_$(input_particle_longname)_input_$(energy_string)keV_$(α)deg_$(number_of_particles)particles.csv /projects/jucl6426/G4EPP/results
    cp /projects/jucl6426/G4EPP/build/results/mlat_65.77deg_input_449.5km_record_450.5km/energy_deposition_$(input_particle_longname)_input_$(energy_string)keV_$(α)deg_$(number_of_particles)particles.csv /projects/jucl6426/G4EPP/results
    """
    )
    close(file)

    global written += 1
  end
end

println("Wrote $(written) files.")