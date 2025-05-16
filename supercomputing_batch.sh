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

# Execute runs
./G4EPP 1000 e- 10.0 0.0
./G4EPP 1000 e- 56.0 0.0
./G4EPP 1000 e- 316.0 0.0
./G4EPP 1000 e- 1778.0 0.0

# Move results out of scratch and rename to correspond to job
cp -r /scratch/alpine/jucl6426/g4epp-build/results/input_450.0km_record_450.0km/* /projects/jucl6426/G4EPP_results/results_$SLURM_JOB_ID
    # Don't need to create the results dir, as the output file argument at the top does that for us

# Exit
echo "Job complete."