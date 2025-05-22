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
  fRad2Deg(3.14159265359 / 180.0),
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
  
  // Set up container for particle properties
  ParticleSample* r = new ParticleSample();
 
  // Start the particle position on a random uniform sampling of a circular area
  G4double theta = G4UniformRand() * 2. * 3.1415926; // u ~ Unif[0, 2 pi)
  G4double radialPosition = G4UniformRand();  // [0, 1)
  G4double diskRadius = 400.*km;

  r->xPos = diskRadius * std::sqrt(radialPosition) * std::cos(theta);
  r->yPos = diskRadius * std::sqrt(radialPosition) * std::sin(theta);
  r->zPos = (fInitialParticleAlt - 500.0)*km; // Subtraction due to coordinate axis location in middle of world volume

  // Particle velocity random variables. Starts electrons with gyromotion about field line at a given pitch angle
  G4double pitchAngle = fBeamPitchAngle * fRad2Deg; // Convert degrees to radians
  G4double gyroPhase  = G4UniformRand() * 2. * fPI;

  // Initial momentum direction of particles in the frame where B points with -z
  double x0 = std::sin(pitchAngle)*std::cos(gyroPhase);
  double y0 = std::sin(pitchAngle)*std::sin(gyroPhase);
  double z0 = -std::cos(pitchAngle);
  
  // Now we need to rotate out of inclined B-field frame into the world frame.
  // We will have B lying in the y-z plane, meaning we need to rotate about the
  // x-axis to get into the world coordinates.
  // Poker Flats is at 65.77 geomagnetic latitude, which has 77.318 deg magnetic 
  // tilt angle => we want to tilt our coordinate system +12.682 deg about x.
  G4double tilt_angle = 12.682 * fRad2Deg; // TODO dip angle will vary based on geomaglat
  r->xDir = x0;
  r->yDir = (std::cos(tilt_angle) * y0) - (std::sin(tilt_angle) * z0);
  r->zDir = (std::sin(tilt_angle) * y0) + (std::cos(tilt_angle) * z0);

  // Verify that pitch angle generation is correct
  double geomagLat_radians = 65.77 * 3.1415926 / 180.0;
  double Bx = 0;
  double By = std::cos(geomagLat_radians);
  double Bz = -2.0 * std::sin(geomagLat_radians);

  double momentumNorm = pow(r->xDir, 2) + pow(r->yDir, 2) + pow(r->zDir, 2);
  double Bnorm = pow(Bx, 2) + pow(By, 2) + pow(Bz, 2);
  double dotProd = (r->xDir * Bx) + (r->yDir * By) + (r->zDir * Bz);
  double generatedPitchAngle_deg = std::acos(dotProd / (momentumNorm * Bnorm)) * 180/3.14159265358979;

  if( abs(generatedPitchAngle_deg - fBeamPitchAngle) > 1){
    G4cout << "** ERROR: Primary generated with pitch angle >1º different than user-specified pitch angle. You should never see this. Please email julia.claxton@colorado.edu with this error and the conditions that produced it." <<
    G4endl << "    generatedPitchAngle_deg: " << generatedPitchAngle_deg <<
    G4endl << "    fBeamPitchAngle: " << fBeamPitchAngle <<
    G4endl;
    throw;
  }

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


