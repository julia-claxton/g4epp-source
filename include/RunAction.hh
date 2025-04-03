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
// $Id: RunAction.hh 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include <vector>
#include "myHistogram.hh"

// Choose your fighter:
// #include "g4root.hh"
//#include "g4xml.hh"
// #include "g4csv.hh"

class G4ParticleDefinition;
class G4Transportation;
class G4CoupledTransportation;
class G4Run;
class SteppingAction;
class myHistogram;

class RunActionMessenger;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void SetHistFileName(G4String name){fHistogramFileName=name;};

    // Helper method to change the Transportation's 'looper' parameters 
    void ChangeLooperParameters(const G4ParticleDefinition* particleDef );

    // Helper method to find the Transportation process for a particle type 
    std::pair<G4Transportation*, G4CoupledTransportation*>
     findTransportation( const G4ParticleDefinition * particleDef,
                         bool reportError= true );

  public:
    void     SetNumberOfTrials( G4int val )   { fNumberOfTrials  = val; }
    void     SetWarningEnergy( double val )   { fWarningEnergy   = val; }
    void     SetImportantEnergy( double val ) { fImportantEnergy = val; }   
    G4int    GetNumberOfTrials()  { return fNumberOfTrials; }
    G4double GetWarningEnergy()   { return fWarningEnergy; }
    G4double GetImportantEnergy() { return fImportantEnergy; } 

  public:
    myHistogram*           fEnergyHist_1;
    myHistogram*           fEnergyHist2D_1;
    myHistogram*           fEnergyHist_2;
    myHistogram*           fEnergyHist2D_2;
  
  private:
    RunActionMessenger*    fRunActionMessenger;
    G4String 		   fHistogramFileName;
  
    // Values for initialising 'loopers' parameters of Transport process
    G4int    fNumberOfTrials  =  0;    // Default will not overwrite
    G4double fWarningEnergy   = -1.0;  // Default values - non operational 
    G4double fImportantEnergy = -1.0;  // Default - will not overwrite

    int    theVerboseLevel = 0;


};


#endif
