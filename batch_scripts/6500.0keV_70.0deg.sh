#!/bin/bash

#SBATCH --job-name G4EPP_6500.0keV_70.0deg
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 23:00:00
#SBATCH --output /projects/jucl6426/G4EPP/results/log_6500.0keV_70.0deg.out
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-20,bhpc-c5-u7-21,bhpc-c5-u7-22,bhpc-c5-u7-23
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Terminate on any non-zero exit status
set -e

# Run simulation
cd /projects/jucl6426/G4EPP/build/
./G4EPP 100000 e- 6500.0 70.0

# Copy results to safe folder
cp -r /projects/jucl6426/G4EPP/build/results/input_450.0km_record_450.0km/backscatter_electron_input_6500.0keV_70.0deg_100000particles.csv /projects/jucl6426/G4EPP/results
cp -r /projects/jucl6426/G4EPP/build/results/input_450.0km_record_450.0km/energy_deposition_electron_input_6500.0keV_70.0deg_100000particles.csv /projects/jucl6426/G4EPP/results

