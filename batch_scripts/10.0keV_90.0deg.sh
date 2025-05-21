#!/bin/bash

#SBATCH --job-name G4EPP_10.0keV_90.0deg
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 05:00:00
#SBATCH --output /projects/jucl6426/G4EPP/results/batch/log_10.0keV_90.0deg.out
#SBATCH --qos=blanca-lair
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Run simulation
cd /projects/jucl6426/G4EPP/build/
./G4EPP 100 e- 10.0 90.0

# Copy results to safe folder
cp -r /projects/jucl6426/G4EPP/build/results/input_450.0km_record_450.0km/backscatter_electron_input_10.0keV_90.0deg_100particles.csv /projects/jucl6426/G4EPP/results/batch
cp -r /projects/jucl6426/G4EPP/build/results/input_450.0km_record_450.0km/energy_deposition_electron_input_10.0keV_90.0deg_100particles.csv /projects/jucl6426/G4EPP/results/batch

