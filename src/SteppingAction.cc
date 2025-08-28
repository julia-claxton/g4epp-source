//
// ********************************************************************
// * License and Disclaimer                       *
// *                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.               *
// *                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.     *
// *                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.      *
// ********************************************************************
//
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingActionMessenger.hh"
#include "G4AutoLock.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"

SteppingAction::SteppingAction(EventAction* eventAction, RunAction* RuAct)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(RuAct),
  fBackscatterFilename(),
  fSteppingMessenger(),
  fCollectionAltitude(450.0)
{
  fSteppingMessenger = new SteppingActionMessenger(this);
}

SteppingAction::~SteppingAction(){delete fSteppingMessenger;}

G4double trackedEnergy = 0;
G4double eEnergy = 0;
std::map<G4int, G4double> allWeights = {{0, 1.0}};

namespace{ G4Mutex aMutex=G4MUTEX_INITIALIZER; } 



void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4AutoLock lock(&aMutex); // Might not be necessary with thread-specific files

  G4Track* track = step->GetTrack();
  G4double trackWeight = track->GetWeight();
  const G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();

  trackedEnergy += step->GetTotalEnergyDeposit() * trackWeight;


  if( particleName == "e-" ){
    eEnergy += step->GetTotalEnergyDeposit() * trackWeight;
  }


  if( step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "eBrem" ){
    G4cout << "Track " << track->GetTrackID() << " produced eBrem. KE = " << step->GetPostStepPoint()->GetKineticEnergy() << G4endl;
  }

  allWeights.insert({track->GetTrackID(), track->GetWeight()});

  //G4double parentWeight = allWeights.at(track->GetParentID());


  if(trackWeight > 0.5){
    G4cout <<
      track->GetTrackID() << "\t" <<
      particleName << "\t" <<
      "Weight: " << trackWeight << "\t" <<
      //"Parent: " << parentWeight << "\t" <<
      track->GetParentID() << "\t" <<
      trackedEnergy / MeV << " MeV" <<
    G4endl;
  }


  
  // ping exits
  if(step->GetTrack()->GetNextVolume() == nullptr){
    G4double exitEnergy = step->GetPostStepPoint()->GetKineticEnergy() * trackWeight;

    if( particleName == "gamma" ){ exitEnergy = step->GetPostStepPoint()->GetTotalEnergy() * trackWeight; }

    G4cout << exitEnergy/keV << " keV " << particleName << " exited world" << G4endl;
    trackedEnergy += exitEnergy;
  }

  // G4cout << eEnergy / MeV << " MeV" << G4endl;
}

void SteppingAction::LogEnergy(G4int histogramAddress, G4double energy){}