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
To build this, run:
mkdir G4EPP-build # Can be anywhere
cd G4EPP-build

cmake -DCMAKE_INSTALL_PREFIX="/Users/luna/Documents/Work/Research/geant4/geant4-install" -DGeant4_DIR="/Users/luna/Documents/Work/Research/geant4/geant4-install/lib" "/Users/luna/Documents/Work/Research/geant4/G4EPP-source"
make

to run:
./G4EPP test_electrons.mac

*/


/// \file electron_detector_main.cc
/// \brief Main program of the electron detector simulation

// Base simulation building classes
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "RunAction.hh"

// Multithreading header support
#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
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

#include <chrono>

// For Printing statistic from Transporation process(es)
#include "G4Electron.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"

// For moving results to subdirectory (because I can't get dataCollection to do it itself
#include <filesystem>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Start simulation timer
  auto t_start = std::chrono::high_resolution_clock::now();

  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

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

  // Construct the default run manager
  #ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    int n_threads = 8; // (Julia's computer)
    runManager->SetNumberOfThreads(n_threads);
    G4cout << "Using " << runManager->GetNumberOfThreads() << " threads." << G4endl;
  #else
    G4RunManager* runManager = new G4RunManager;
  #endif

  // Physics list
  G4PhysListFactory factory;
  G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QBBC");
  
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  physicsList->SetVerboseLevel(0);
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization());

  G4double lowLimit = 250. * eV;
  G4double highLimit = 100. * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);


  // Fine grained control of thresholds for looping particles
  auto runAction= new RunAction();
  runAction->SetWarningEnergy(   1.0 * keV );
              // Looping particles with E < ^ keV will be killed after 1 step
              //   with warning.
              // Looping particles with E > ^ keV will generate a warning.
  runAction->SetImportantEnergy( 0.1 * keV );
  runAction->SetNumberOfTrials( 50 ); // Looping particles with E > 0.1 keV will survive for up to 30 'tracking' steps, and only be killed if they still loop.
  // Note: this mechanism overwrites the thresholds established by the call to UseLowLooperThresholds() above.
  
  runManager->SetUserAction(runAction);


  // Suppress large verbosity from EM & hadronic processes
  G4EmParameters::Instance()->SetVerbose(-1);
  G4HadronicProcessStore::Instance()->SetVerbose(0);


  // Initialize visualization
  /*
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  visManager->Initialize();
  */

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // Batch mode
    G4String execute = "/control/execute ";
    G4String command_to_run;

    // If one argument is provided, interpret it as a macro to run
    if(argc == 2){
      command_to_run = argv[1];
      std::cout << "Running in macro mode with " + command_to_run << std::endl;
    }

    // If 2 arguments are provided, interpret them as an energy and pitch angle to run
    if(argc == 3){
      G4String energy = argv[1];
      G4String pitch_angle = argv[2];
      UImanager->ApplyCommand("/control/alias BEAM_ENERGY_KEV " + energy);
      UImanager->ApplyCommand("/control/alias BEAM_PITCH_ANGLE " + pitch_angle);

      std::cout << "Running in single beam mode at " + energy + " keV, " + pitch_angle + " deg" << std::endl;
      command_to_run = "run_single_beam.mac";
    }

    UImanager->ApplyCommand(execute + command_to_run);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute initialize_visualization.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  // delete visManager;
  delete runManager;

  /*
  // Move output files to results directory. Can't get dataCollection to do this automatically so we resort to system commands.
  std::cout << "Moving output files to ./results..." << std::endl;
  system("mv -f $(find . -name \"electron_*\") ./results");
  system("mv -f $(find . -name \"photon_*\") ./results");
  */

  // End simulation timer
  auto t_end = std::chrono::high_resolution_clock::now();

  // Calculate elapsed simulation (wall) time
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  std::cout << "Simulation completed in : " << elapsed_time_ms/1000. << " seconds" << std::endl;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....