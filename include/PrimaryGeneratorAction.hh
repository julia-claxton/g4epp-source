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
// $Id: PrimaryGeneratorAction.hh 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"

class G4ParticleGun;
// class G4GeneralParticleSource;
class G4Event;
class G4Box;
class PrimaryGeneratorMessenger;

struct ParticleSample{
	G4double xPos, yPos, zPos;
	G4double xDir, yDir, zDir;
	G4double energy;
};

/// The primary generator action class with particle gun
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         

    void GenerateElectrons(ParticleSample*);
    void GenerateSolarSpectra(ParticleSample*);
    void GenerateCXB(ParticleSample*);
    
    // Messenger methods
    void SetBeamEnergy(G4double E0){ fBeamEnergy = E0;};
    void SetBeamPitchAngle(G4double maxPA){fBeamPitchAngle = maxPA; };
    void SetParticleInitialAlt(G4double alt){fInitialParticleAlt = alt; };
    void SetInputParticleType(G4String type){fSourceType = type; };
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; } // method to access particle gun
  
  private:
    G4ParticleGun*  fParticleGun; // Pointer to G4 gun class
    PrimaryGeneratorMessenger* fPrimaryMessenger;
    G4double fBeamEnergy;
    G4double fBeamPitchAngle;
    G4double fInitialParticleAlt;
    G4double fPI;
    G4double fRad2Deg;
    G4String fSourceType;
};

#endif