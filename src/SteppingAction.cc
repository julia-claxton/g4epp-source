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
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
// #include "DetectorAnalysis.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingActionMessenger.hh"
#include "G4AutoLock.hh"

// Initialize autolock for multiple threads writing into a single file
namespace{ G4Mutex aMutex=G4MUTEX_INITIALIZER; } 

SteppingAction::SteppingAction(EventAction* eventAction, RunAction* RuAct)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(RuAct),
  fEnergyThreshold_keV(0.),
  //fWindowAlt(500.),
  fDataCollectionType(0),
  fBackscatterFilename(),
  fSteppingMessenger()
{
  fSteppingMessenger = new SteppingActionMessenger(this);
}

SteppingAction::~SteppingAction(){delete fSteppingMessenger;}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // ===========================
  // Guard Block
  // ===========================
  G4Track* track = step->GetTrack();

  // Check for NaN energy
  if( std::isnan(step->GetPostStepPoint()->GetKineticEnergy()) )
  {  
    track->SetTrackStatus(fStopAndKill);
    G4cout << "Particle killed for negative energy." << G4endl;
    G4cout << "Particle killed at: " << step->GetPreStepPoint()->GetKineticEnergy()/keV << " keV , Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
  }

  // Check for exceeding 1 second of simulation time
  G4double time = track->GetProperTime();
  if(time/second > 1)
  {
    track->SetTrackStatus(fStopAndKill);
    G4cout << "Particle killed for time > 1s." << G4endl;
    G4cout << "Particle killed at: " << step->GetPreStepPoint()->GetKineticEnergy()/keV << " keV , Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
  }

  // ===========================
  // Ionization Tracking
  // ===========================

  // TODO

  // TODO headers on result files


  // ===========================
  // Backscatter Tracking
  // ===========================
  const G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
  const G4ThreeVector momentum = track->GetMomentumDirection();
  const G4ThreeVector position = track->GetPosition();
  const G4double energy =  step->GetPreStepPoint()->GetKineticEnergy();
        
  if( (position.z()/km > 450.0-500.0) && (momentum.z() > 0) ) // If particle is above 450 km and moving upwards. Subtract 500 because 0.0 in world coordinates = +500 km above sea level
  {
    // Lock scope to stop threads from overwriting data in same file
    G4AutoLock lock(&aMutex);

    // Write position and momentum to file
    std::ofstream dataFile;
    dataFile.open(fBackscatterFilename, std::ios_base::app); // Open file in append mode
    dataFile 
      << particleName << ','
      << position.x()/m << ',' 
      << position.y()/m << ','
      << (position.z()/m) + 500000 << ',' // Shift so we are writing altitude above sea level to file rather than the world coordinates
      << momentum.x() * energy/keV << ','
      << momentum.y() * energy/keV << ','
      << momentum.z() * energy/keV << '\n'; 
    dataFile.close();

    // Kill particle after data collection
    track->SetTrackStatus(fStopAndKill);
    G4cout << "Recorded and killed upgoing " << particleName << " at 450 km." << G4endl; // Status message
  }



  /*
  switch(fDataCollectionType)
  {
    
    case(0):  // Collects energy deposition per altitude
    { 
    	// Gets energy delta of particle over step length
    	const G4double energyDep = step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy();
          
      if(energyDep > fEnergyThreshold_keV*keV)
      {
        // Gets altitude of particle
        G4ThreeVector position = track->GetPosition();
        G4double      zPos     = position.z();

        // Adds energy deposition to vector owned by RunAction, which is
        // written to a results file per simulation run
        G4int altitudeAddress = std::floor(500. + zPos/km);

        // Thread lock this so only one thread can deposit energy into
        // the histogram at a time. Unlocks when l goes out of scope.
        if(altitudeAddress > 0 && altitudeAddress < 1000) 
        {
          LogEnergy(altitudeAddress, energyDep/keV);
        }
      }

      // If not e- or gamma, move on
      break;
    }

    case(4): // Radiation and ionization histograms 
    {
      // Electron analysis
      if(track->GetDynamicParticle()->GetDefinition()->GetParticleName() == "e-")
      {
        // Gets energy delta of particle over step length
        G4double energyBefore = step->GetPreStepPoint()->GetKineticEnergy(); 
        G4double energyAfter = step->GetPostStepPoint()->GetKineticEnergy();
        G4double energyDep = energyBefore - energyAfter;

        // Gets altitude of particle
        G4ThreeVector position = track->GetPosition();
        G4double      zPos     = position.z();
    
        // Adds energy deposition to vector owned by RunAction, which is
        // written to a results file per simulation run
        G4int altitudeAddress = std::floor(500. + zPos/km);
        
        // Check for valid altitude address
        if(altitudeAddress > 0 && altitudeAddress < 1000 && energyDep > 0. && energyAfter > 0. && !std::isnan(energyDep) && !std::isnan(energyAfter)) 
          // check if energies are nan 
        {
          LogEnergyToSpecificHistogram(altitudeAddress, energyDep, energyAfter, 1);
        }
      }

      // Check if particle is a photon 
      else if(track->GetDynamicParticle()->GetDefinition()->GetParticleName() == "gamma")
      {
        // Gets energy delta of particle over step length
        G4double energyBefore = step->GetPreStepPoint()->GetKineticEnergy(); 
        G4double energyAfter = step->GetPostStepPoint()->GetKineticEnergy();
        G4double energyDep = energyBefore - energyAfter;
          
        // Gets altitude of particle
              G4ThreeVector position = track->GetPosition();
              G4double      zPos     = position.z();
            
              // Adds energy deposition to vector owned by RunAction, which is
              // written to a results file per simulation run
        G4int altitudeAddress = std::floor(500. + zPos/km);
      
        // Check for valid altitude address
        if(altitudeAddress > 0 && altitudeAddress < 1000 && energyDep > 0. && energyAfter > 0. && !std::isnan(energyDep) && !std::isnan(energyAfter)) 
        {
          LogEnergyToSpecificHistogram(altitudeAddress, energyDep, energyAfter, 2);
        }
      }
      
      // If not e- or gamma, move on
      break;
    }
    
    default: 
      throw std::runtime_error("Enter a valid data collection type!");
      break;
  }
  */
}

void SteppingAction::LogEnergy(G4int histogramAddress, G4double energy)
{
  G4AutoLock lock(&aMutex);
  fRunAction->fEnergyHist_1->AddCountToBin(histogramAddress, energy/keV);
}

void SteppingAction::LogEnergyToSpecificHistogram(G4int histogramAddress, G4double entry1, G4double entry2, G4int whichHistogram)
{
  G4AutoLock lock(&aMutex);
  switch(whichHistogram)
  {
    case(1):
      fRunAction->fEnergyHist_1->AddCountToBin(histogramAddress, entry1/keV);
      fRunAction->fEnergyHist2D_1->AddCountTo2DHistogram(histogramAddress, entry2/keV);
      break;
    case(2):
      fRunAction->fEnergyHist_2->AddCountToBin(histogramAddress, entry1/keV);
      fRunAction->fEnergyHist2D_2->AddCountTo2DHistogram(histogramAddress, entry2/keV);
      break;
    default:
      throw std::runtime_error("Enter a valid histogram selection!");
      break;
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
