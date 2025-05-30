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
// $Id: SteppingAction.hh 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>


class EventAction;
class RunAction;
class G4LogicalVolume;

class SteppingActionMessenger;

/// Stepping action class
///

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(EventAction* eventAction, RunAction* RuAct);
    virtual ~SteppingAction();

    void SetBackscatterFilename(G4String name){fBackscatterFilename = name;};
    void LogEnergy(G4int, G4double);    
    
    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

    // Messenger methods
    void SetCollectionAltitude(G4double collectionAltitude){ fCollectionAltitude = collectionAltitude;};

  private:
    EventAction*             fEventAction;
    RunAction*               fRunAction;
    G4double                 fEnergyThreshold_keV;
    G4String                 fBackscatterFilename;
    SteppingActionMessenger* fSteppingMessenger;
    G4double                 fCollectionAltitude;
};


#endif
