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
// $Id: RunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"
#include "G4Electron.hh"

#include "myHistogram.hh"
#include "RunActionMessenger.hh"

#include <fstream>
#include "../include/csv.h"

RunAction::RunAction():
  G4UserRunAction(),
  fRunActionMessenger(),
  fEnergyDepositionFileName(),
  fBackscatterFilename()
{
  fWarningEnergy = 0.01 * keV; // Particles below this energy are killed after 1 step. Value arbitrary 
  fImportantEnergy = 0.1 * keV; // Particles above this energy are killed after fNumberOfTrials if they are looping. Value arbitrary 
  fNumberOfTrials = 1000; // Number of trials before a looping 'important' particle is killed. Value arbitrary. Set very high to avoid losing particles

  fRunActionMessenger = new RunActionMessenger(this); 

  // Initialize energy deposition histogram
  fEnergyDepositionHistogram = new myHistogram(); // Defaults to 0-999 km in 1 km bins
}

RunAction::~RunAction()
{
  delete fEnergyDepositionHistogram;
  delete fRunActionMessenger;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  int threadID = G4Threading::G4GetThreadId();

  // If we are the main thread, create backscatter file and write header
  if(threadID == -1)
  {
    // after adding this block, writes in steppingaction stopped working. fix!!
    std::ofstream dataFile;
    dataFile.open(fBackscatterFilename, std::ios_base::app); // Open file in append mode
    dataFile << "particle_name,x_meters,y_meters,z_meters,px,py,pz\n"; // TODO momentum units??
    dataFile.close();
    G4cout << "datafile closed" << G4endl;
  }

  // Change parameters of looping particles. TODO delete?
  ChangeLooperParameters( G4Electron::Definition() );
}

void RunAction::ChangeLooperParameters(const G4ParticleDefinition* particleDef )
{
  if( particleDef == nullptr )
    particleDef = G4Electron::Definition();
  auto transportPair= findTransportation( particleDef );
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;

  if( transport != nullptr )
  { 
    // Change the values of the looping particle parameters of Transportation 
    if( fWarningEnergy   >= 0.0 )
      transport->SetThresholdWarningEnergy(  fWarningEnergy ); 
    if( fImportantEnergy >= 0.0 )
      transport->SetThresholdImportantEnergy(  fImportantEnergy ); 
    if( fNumberOfTrials  > 0 )
      transport->SetThresholdTrials( fNumberOfTrials );
  }
  else if( coupledTransport != nullptr )
  { 
    // Change the values for Coupled Transport
    if( fWarningEnergy   >= 0.0 )
      coupledTransport->SetThresholdWarningEnergy(fWarningEnergy); 
    if( fImportantEnergy >= 0.0 )
      coupledTransport->SetThresholdImportantEnergy(fImportantEnergy); 
    if( fNumberOfTrials  > 0 )
      coupledTransport->SetThresholdTrials(fNumberOfTrials);
  }
}


void RunAction::EndOfRunAction(const G4Run*)
{
  int threadID = G4Threading::G4GetThreadId();

  // If we are not the main thread, write energy deposition to file and exit
  if(threadID != -1)
  {
    std::string threadFilename = fEnergyDepositionFileName.substr(0, fEnergyDepositionFileName.length()-4) + "_thread" + std::to_string(threadID) + ".csv"; // Thread-specific filename
    fEnergyDepositionHistogram->WriteHistogramToFile(threadFilename);
    G4cout << "Thread " + std::to_string(threadID) + " complete" << G4endl;
    return;
  }

  // If we are the main thread, merge energy deposition datafiles from each thread. Main thread ends after workers are done, so this is the end of the simulation
  G4cout << "Main Thread: Merging energy deposition data... ";

  // Instantiate merged histogram
  myHistogram* mainEnergyDepositionHistogram = new myHistogram(); // 1000 km in 1 km bins

  // Add energy deposition from each thread to the merged histogram
  for(int threadFileToMerge = 0; threadFileToMerge < 8; threadFileToMerge++)
  {
    std::string threadFilename = fEnergyDepositionFileName.substr(0, fEnergyDepositionFileName.length()-4) + "_thread" + std::to_string(threadFileToMerge) + ".csv"; // Thread-specific filename

    // Read in CSV from this thread
    io::CSVReader<2> in(threadFilename);
    in.read_header(io::ignore_extra_column, "altitude_km", "energy_deposition_kev");
    int altitudeAddress; double energy_deposition;

    // For each row in the file, add energy deposition to main histogram
    while(in.read_row(altitudeAddress, energy_deposition)){ 
      mainEnergyDepositionHistogram->AddCountToBin(altitudeAddress, energy_deposition); // TODO isn't working in for loop -- why?
    }

    // Delete this thread-specific file
    std::remove(threadFilename.c_str());
  }

  // Write main energy histogram to file
  mainEnergyDepositionHistogram->WriteHistogramToFile(fEnergyDepositionFileName);

  // Status message and exit
  G4cout << G4endl << "Done" << G4endl;
}


std::pair<G4Transportation*, G4CoupledTransportation*> RunAction::findTransportation(const G4ParticleDefinition* particleDef, bool reportError)
{
  const auto *partPM=  particleDef->GetProcessManager();
    
  G4VProcess* partTransport = partPM->GetProcess("Transportation");
  auto transport= dynamic_cast<G4Transportation*>(partTransport);

  partTransport = partPM->GetProcess("CoupledTransportation");
  auto coupledTransport=
     dynamic_cast<G4CoupledTransportation*>(partTransport);

  if( reportError && !transport && !coupledTransport )
  {
     G4cerr << "Unable to find Transportation process for particle type "
            << particleDef->GetParticleName()
            << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
            << G4endl;
  }
  
  return std::make_pair( transport, coupledTransport );
}


