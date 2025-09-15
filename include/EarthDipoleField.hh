#ifndef EARTHDIPOLEFIELD_h
#define EARTHDIPOLEFIELD_h 1

#include "G4MagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class EarthDipoleFieldMessenger;

class EarthDipoleField : public G4MagneticField
{
public:
  EarthDipoleField();
  virtual ~EarthDipoleField() override;

  // From base class
  virtual void GetFieldValue(const G4double Point[4], G4double *Bfield) const override;

  // Messenger methods
  void SetMLAT(G4double MLAT_deg){ fMLAT_degrees = MLAT_deg; };

private:
  EarthDipoleFieldMessenger* fDipoleFieldMessenger;
  G4double fDipoleMoment;
  G4double fMLAT_degrees;
  G4double fRe;
  G4double fu0;
  G4double fpi;
  G4double fEarthRadius; 
};

#endif
