#!/bin/bash

#SBATCH --job-name G4EPP
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 50:00:00
#SBATCH --output /projects/jucl6426/G4EPP_results/results_%j/log.out
#SBATCH --qos=blanca-lair
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Print slurm parameters
echo "
#SBATCH --job-name G4EPP
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 50:00:00
#SBATCH --output /projects/jucl6426/G4EPP_results/results_%j/log.out
#SBATCH --qos=blanca-lair
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu
"

# https://stackoverflow.com/questions/65603381/slurm-nodes-tasks-cores-and-cpus
# https://slurm.schedmd.com/sbatch.html

# Remove previous results
rm -f /scratch/alpine/jucl6426/g4epp-build/results/*
cd /scratch/alpine/jucl6426/g4epp-build

# Execute runs
for e in 10.0 56.0 316.0 1778.0 10000.0; do
  for pa in 0.0 18.0 36.0 54.0 72.0 90.0; do
    # Run simulation
    ./G4EPP 100 e- $e $pa
    # Move results
    cp -r /projects/jucl6426/G4EPP/build/results/input_450.0km_record_450.0km/* /projects/jucl6426/G4EPP/results/results_$SLURM_JOB_ID
  done
done

# Exit
echo "Job complete."