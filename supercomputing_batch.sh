#!/bin/bash

#SBATCH --job-name G4EPP
#SBATCH --nodes 1
#SBATCH --threads-per-core 1
#SBATCH --ntasks-per-node 128
#SBATCH --time 00:5:00
#SBATCH --output /projects/jucl6426/G4EPP_results/G4EPP_LOG_%j.txt
#SBATCH --qos=preemptable
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Remove previous results
rm -f /scratch/alpine/jucl6426/g4epp-build/results/*

# Execute runs
./G4EPP 1000 e- 10000 67

# Move results out of scratch
rm -r /projects/jucl6426/G4EPP_output/* # Clear old results. Will fail if doesn't exist
cp -r /scratch/alpine/jucl6426/g4epp-build/results /projects/jucl6426/G4EPP_results

# Exit
echo "Job complete."