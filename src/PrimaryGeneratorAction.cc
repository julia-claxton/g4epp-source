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
// $Id: PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "PrimaryGeneratorMessenger.hh"
#include <math.h>


PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fPrimaryMessenger(0),
  fEnergyDistType(0),
  fPitchAngleDistType(0),
  fE0(100.),
  fMaxPitchAngle(40.),
  fInitialParticleAlt(500.),
  fPI(3.14159265359),
  fRad2Deg(180. / 3.14159265359),
  fSourceType(0)
{
  fParticleGun  = new G4ParticleGun();
  
  fPrimaryMessenger = new PrimaryGeneratorMessenger(this);
  
  G4ParticleDefinition* electronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(electronParticle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateElectrons(ParticleSample* r)
{

  G4String word;
  std::ifstream inputFile; 

  // Initial position RV's
  G4double theta = G4UniformRand() * 2. * 3.1415926; // u ~ Unif[0, 2 pi)
  G4double radialPosition = G4UniformRand();  // [0, 1)
  G4double diskRadius = 400.*km;

  // Random uniform sampling on a circular area
  r->xPos = diskRadius * std::sqrt(radialPosition) * std::cos(theta);
  r->yPos = diskRadius * std::sqrt(radialPosition) * std::sin(theta);
  
  // Subtraction due to coordinate axis location in middle of volume
  std::cout << "Initial altitude: " << fInitialParticleAlt << " km" << std::endl;
  r->zPos = (fInitialParticleAlt - 500)*km;

  // Particle attribute random variables
  // Starts electrons with gyromotion about field line
  G4double maxPitchAngle = fMaxPitchAngle * 3.1415926 / 180.;   // rad
  
  // Angular RV's
  G4double gyroPhase  = G4UniformRand() * 2. * 3.1415926;
  G4double pitchAngle = maxPitchAngle;

  // NB: need to rotate into inclined B-field frame
  // Poker Flats: 65.77 geomagnetic latitude --> 77.318 deg magnetic tilt angle
  // => we want to tilt our coordinate system 12.682 deg in the z-y plane
  G4double tilt_angle = 12.682 * fPI / 180.; // rad

  // Initial momentum direction of particles
  r->xDir = std::sin(pitchAngle)*std::cos(gyroPhase);
  r->yDir = std::sin(pitchAngle)*std::sin(gyroPhase);
  r->zDir = -std::cos(pitchAngle);

  // Rotate
  r->yDir = std::cos(tilt_angle) * r->yDir - std::sin(tilt_angle) * r->zDir;
  r->zDir = std::sin(tilt_angle) * r->yDir + std::cos(tilt_angle) * r->zDir;

  // Set energy
  r->energy = fE0 * keV;
}

void PrimaryGeneratorAction::GenerateSolarSpectra(ParticleSample* r)
{

}

void PrimaryGeneratorAction::GenerateCXB(ParticleSample* r)
{

}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Struct to store particle initial properties: 
  // position           - (x, y, z) [km]
  // momentum direction - (px, py, pz) []
  // energy 		- E [keV]
  ParticleSample* r = new ParticleSample();
 
  // Generator method to fill ParticleSample struct
  switch(fSourceType)
  {
    case(0):
      GenerateElectrons(r);
      break;
 
    case(1):
      GenerateSolarSpectra(r);
      break;
    
    case(2):
      GenerateCXB(r);
      break;
    
    default:
      throw std::invalid_argument("Enter a valid source type");
  }

  fParticleGun->SetParticlePosition(
		  G4ThreeVector(r->xPos, r->yPos, r->zPos)); 
  
  fParticleGun->SetParticleMomentumDirection(
		  G4ThreeVector(r->xDir, r->yDir, r->zDir));
  
  fParticleGun->SetParticleEnergy(r->energy);
  
  // Geant method to create initial particle with the above properties 
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // Free memory from ParticleSample struct
  delete r;
}


