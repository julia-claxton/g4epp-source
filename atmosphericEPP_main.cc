//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

/*
# To build the G4EPP executable: (in zsh)
cd path/to/G4EPP-build # You'll need to create this directory if it doesn't exist already
rm -r /path/to/G4EPP-build/ *  # Remove space between / and * at the end before running (C++ gets angry if there's a slash-star in a block comment). Removes any existing files or old builds
cmake -DCMAKE_INSTALL_PREFIX="/path/to/geant4-install" -DGeant4_DIR="/path/to/geant4-install/lib" /path/to/G4EPP-source
make # Build executable
chmod +x ./RUN_ALL.sh # Make the runall script executable by anyone

# To run the simulation
./G4EPP <BEAM ENERGY IN KEV> <BEAM PITCH ANGLE IN DEG>
# or
./RUN_ALL.sh
*/

/// \file atmosphericEPP_main.cc
/// \brief Main function to run Geant4 EPP simulation

// Base simulation building classes
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "RunAction.hh"

// Multithreading support
#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
  #include "G4Threading.hh"
#else
  #include "G4RunManager.hh"
#endif

// Physics lists
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4EmParameters.hh"
#include "G4HadronicProcessStore.hh"

// For time display
#include <chrono>

// For Printing statistic from Transporation process(es)
#include "G4Electron.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"

int main(int argc,char** argv)
{
  // Start simulation timer
  auto t_start = std::chrono::high_resolution_clock::now();

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Set the seeds
  long seeds[2];
  time_t systime = time(NULL);
  
  // Seed built in c-rand engine
  srand (systime);

  // Geant rand engine
  seeds[0] = (long) systime;
  seeds[1] = (long) (systime*G4UniformRand());
  G4Random::setTheSeeds(seeds);

  // Construct the run manager
  #ifdef G4MULTITHREADED
    // Multithreaded mode
    G4MTRunManager* runManager = new G4MTRunManager;
    int n_threads = 40; // Change this number to the desired number of threads
    runManager->SetNumberOfThreads(n_threads);
    G4cout << "Using " << runManager->GetNumberOfThreads() << " threads." << G4endl;
  #else
    // Singlethreaded mode
    G4RunManager* runManager = new G4RunManager;
  #endif

  // Physics list
  G4PhysListFactory factory;
  G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QBBC"); // QBBC uses EM v1. Do we need more updated EM model? TODO
  
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  physicsList->SetVerboseLevel(0);
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization());

  G4double lowLimit = 250. * eV;
  G4double highLimit = 100. * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);

  // Suppress large verbosity from EM & hadronic processes
  G4EmParameters::Instance()->SetVerbose(-1);
  G4HadronicProcessStore::Instance()->SetVerbose(0);

  // Get the pointer to the user interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Decide what to execute based on command line input
  G4String execute = "/control/execute ";
  G4String command_to_run;
  
  // If one argument is provided, interpret it as a macro to run
  if(argc == 2){
    command_to_run = argv[1];
    std::cout << "=====================================================================" << std::endl;
    std::cout << "Running in macro mode with " + command_to_run << std::endl;
    std::cout << "=====================================================================" << std::endl;
  }
  // If 3 arguments are provided, interpret them as a particle, energy, and pitch angle to run
  else if(argc == 4){
    // Set particle definition variable (uses Geant4's particle names: https://fismed.ciemat.es/GAMOS/GAMOS_doc/GAMOS.5.1.0/x11519.html
    G4String particle = argv[1];
    UImanager->ApplyCommand("/control/alias BEAM_PARTICLE " + particle); 
    
    // Set particle longname - what the result file will call the input particle. This is just for clarity to 
    // the end user on what each result file represents, as I think "photon" and "electron" are clearer than
    // G4's internal names "gamma" and "e-".
    G4String longname;
    if(particle == "e-")         {longname = "electron";}
    else if(particle == "gamma") {longname = "photon";}
    else                         {longname = particle;}
    UImanager->ApplyCommand("/control/alias BEAM_PARTICLE_LONGNAME " + longname);

    // Set energy and pitch angle of beam
    G4String energy = argv[2];
    G4String pitch_angle = argv[3];
    UImanager->ApplyCommand("/control/alias BEAM_ENERGY_KEV " + energy);
    UImanager->ApplyCommand("/control/alias BEAM_PITCH_ANGLE_DEG " + pitch_angle);

    std::cout << "=====================================================================" << std::endl;
    std::cout << "Running in single beam mode at " + energy + " keV, " + pitch_angle + " deg" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    command_to_run = "run_single_beam.mac";
  }
  else
  {
    std::cout << "Incorrect number of command line arguments provided!" << std::endl;
    throw;
  }

  // Report current system time before starting the run
  std::time_t current_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  G4cout << "Starting simulation: " << std::ctime(&current_time) << G4endl;
  
  // Execute run
  UImanager->ApplyCommand(execute + command_to_run);

  // End run
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  delete runManager;

  // End simulation timer
  auto t_end = std::chrono::high_resolution_clock::now();

  // Report elapsed simulation time (realtime, not simulation time)
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  std::cout << "Simulation completed in : " << elapsed_time_ms/1000.0 << " seconds" << std::endl;

  return 0;
}