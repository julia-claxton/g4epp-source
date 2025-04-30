#!/bin/bash

#SBATCH --job-name G4EPP
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 24:00:00
#SBATCH --output /projects/jucl6426/G4EPP_results/G4EPP_LOG_%j.txt
#SBATCH --qos=blanca-lair
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Remove previous results
rm /scratch/alpine/jucl6426/g4epp-build/results/*

# Execute runs
./G4EPP e- 10000 67

# Move results out of scratch
rm -r /projects/jucl6426/G4EPP_output/* # Clear old results. Will fail if doesn't exist
cp -r /scratch/alpine/jucl6426/g4epp-build/results /projects/jucl6426/G4EPP_results

# Exit
echo "Job complete."