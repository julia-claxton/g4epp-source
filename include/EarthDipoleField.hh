

#ifndef EARTHDIPOLEFIELD_h
#define EARTHDIPOLEFIELD_h 1

#include "G4MagneticField.hh"
#include "G4SystemOfUnits.hh"

class EarthDipoleField : public G4MagneticField
{
public:

  EarthDipoleField();
  virtual ~EarthDipoleField() override;

  virtual void GetFieldValue(const G4double Point[4], 
		                   G4double *Bfield) const override;
private:
  G4double fDipoleMoment;
  G4double fGeomagLatitude;
  G4double fEarthRadius; 

};

#endif
