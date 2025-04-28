#!/bin/bash

#SBATCH --job-name G4EPP
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 24:00:00
#SBATCH --output G4EPP_LOG_%j.txt
#SBATCH --qos=blanca-lair
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

###############################################
#
# START FROM HERE IF BUILDING G4EPP FROM SOURCE
#
###############################################
cd /projects/$USER
cd g4epp-source
git pull

# Create/reset build directory
cd /scratch/alpine/$USER
rm -r /scratch/alpine/$USER/g4epp-build # May fail if already doesn't exist yet
mkdir /scratch/alpine/$USER/g4epp-build

# Build executable
cd g4epp-build
cmake -DCMAKE_INSTALL_PREFIX="/projects/$USER/geant4/install" -DGeant4_DIR="/projects/$USER/geant4/install/lib64" "/projects/$USER/g4epp-source"
make

# Make batch script executable
chmod +x ./supercomputing_batch.sh

# Remove previous results
rm /scratch/alpine/jucl6426/g4epp-build/results/*

# Execute runs
./G4EPP e- 10000 67

# Move results out of scratch
rm -r /projects/jucl6426/G4EPP_output/* # Clear old results. Will fail if doesn't exist
cp -r /scratch/alpine/jucl6426/g4epp-build/results /projects/jucl6426/G4EPP_results

# Exit
echo "Job complete."