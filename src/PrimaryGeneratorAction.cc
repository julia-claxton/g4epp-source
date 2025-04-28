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
/// \brief Creates particles that will be propagated in the simulation

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
  fBeamEnergy(100.),
  fBeamPitchAngle(40.),
  fInitialParticleAlt(450.0),
  fPI(3.14159265359),
  fRad2Deg(180. / 3.14159265359),
  fSourceType("e-")
{
  fPrimaryMessenger = new PrimaryGeneratorMessenger(this);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Select input particle type
  fParticleGun  = new G4ParticleGun();
  G4ParticleDefinition* inputParticle = G4ParticleTable::GetParticleTable()->FindParticle(fSourceType); // Electron = "e-", proton = "proton", photon = "gamma"
  fParticleGun->SetParticleDefinition(inputParticle);
  
  // Struct to store particle initial properties: 
  // position           - (x, y, z) [km]
  // momentum direction - (px, py, pz) []
  // energy 		- E [keV]
  ParticleSample* r = new ParticleSample();
 
  // Initial position random variables
  // Starts the particle position on a random uniform sampling of a circular area
  G4double theta = G4UniformRand() * 2. * 3.1415926; // u ~ Unif[0, 2 pi)
  G4double radialPosition = G4UniformRand();  // [0, 1)
  G4double diskRadius = 400.*km;

  r->xPos = diskRadius * std::sqrt(radialPosition) * std::cos(theta);
  r->yPos = diskRadius * std::sqrt(radialPosition) * std::sin(theta);
  r->zPos = (fInitialParticleAlt - 500.0)*km; // Subtraction due to coordinate axis location in middle of world volume

  // Particle velocity random variables. Starts electrons with gyromotion about field line
  G4double pitchAngle = fBeamPitchAngle * fPI / 180.0; // Convert degrees to radians
  G4double gyroPhase  = G4UniformRand() * 2. * 3.1415926;

  // Initial momentum direction of particles
  // TODO collection altitude higher?
  r->xDir = std::sin(pitchAngle)*std::cos(gyroPhase);
  r->yDir = std::sin(pitchAngle)*std::sin(gyroPhase);
  r->zDir = -std::cos(pitchAngle);

  // Need to rotate into inclined B-field frame since B-field is inclined from the world z axis.
  // Poker Flats: 65.77 geomagnetic latitude has 77.318 deg magnetic tilt angle => we want to 
  // tilt our coordinate system 12.682 deg in the z-y plane.
  G4double tilt_angle = 12.682 * fPI / 180.0; // Convert degrees to radians
  r->yDir = std::cos(tilt_angle) * r->yDir - std::sin(tilt_angle) * r->zDir;
  r->zDir = std::sin(tilt_angle) * r->yDir + std::cos(tilt_angle) * r->zDir;

  // Set energy
  r->energy = fBeamEnergy * keV;

  // Communicate to particle gun
  fParticleGun->SetParticlePosition(G4ThreeVector(r->xPos, r->yPos, r->zPos)); 
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(r->xDir, r->yDir, r->zDir));
  fParticleGun->SetParticleEnergy(r->energy);
  
  // Geant method to create initial particle with the above properties 
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // Free memory from ParticleSample struct
  delete r;
}


