#!/bin/bash

#SBATCH --job-name G4EPP_238.0keV_75.0deg
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 7-00:00:00
#SBATCH --output /projects/jucl6426/G4EPP/results/log_238.0keV_75.0deg.out
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-22,bhpc-c5-u7-23
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Terminate on any non-zero exit status
set -e

# Run simulation
cd /projects/jucl6426/G4EPP/build/
./G4EPP 100000 e- 238.0 75.0

# Copy results to safe folder
cp /projects/jucl6426/G4EPP/build/results/input_449.5km_record_450.5km/backscatter_electron_input_238.0keV_75.0deg_100000particles.csv /projects/jucl6426/G4EPP/results
cp /projects/jucl6426/G4EPP/build/results/input_449.5km_record_450.5km/energy_deposition_electron_input_238.0keV_75.0deg_100000particles.csv /projects/jucl6426/G4EPP/results

