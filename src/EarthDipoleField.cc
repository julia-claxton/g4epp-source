

#include "EarthDipoleField.hh"

/* The following class calculates the Earth's magnetic field strength 
 * and direction according to a tilted dipole model. This class inherits
 * from G4MagneticField. GetFieldValue() is a Geant virtual method that 
 * is called to obtain magnetic and electric field values for particle 
 * propagation purposes.
 * 
 * Future work: 
 *  - implement IGRF magnetic field (and others)
 * 
 * Grant Berland
 */


EarthDipoleField::EarthDipoleField()
: G4MagneticField(),
  fDipoleMoment(8.05e6),    // Tesla-km^3
  fGeomagLatitude(65.77),   // deg, Poker Flat geomagnetic latitude
  fEarthRadius(6371.)       // km
{}


EarthDipoleField::~EarthDipoleField()
{}


void EarthDipoleField::GetFieldValue(const G4double Point[4],
		                             G4double *Bfield) const
{
  // Point is a spacetime 4-vector:
  // Point[0..3] ~ (x, y, z, t)
  // Bfield is a pointer to a 6x1 array of E- and B-field components 
	
  G4double geomagLat_radians = fGeomagLatitude * 3.1415926 / 180.;

  // Radial distance from Earth center, input in kilometers
  // 1000 km / 2 = 500 km addition to account for coordinate axes in 
  // center of simulation volume
  G4double z = fEarthRadius + (Point[2]/km + 1000./2.);  // km

  // Magnitude of B-field, units assigned here
  G4double B_magnitude = fDipoleMoment / std::pow(z, 3) * tesla; // T

  // Bfield[0] ~ West direction
  // Bfield[1] ~ North direction
  // Bfield[2] ~ Up direction, or radially out from Earth 
  Bfield[0] = 0; // Bx
  Bfield[1] = B_magnitude * std::cos(geomagLat_radians);       // By
  Bfield[2] = B_magnitude * -2. * std::sin(geomagLat_radians); // Bz
  Bfield[3] = 0; // Ex
  Bfield[4] = 0; // Ey
  Bfield[5] = 0; // Ez

  // Debugging print statement
  //std::cout << "Bz = " << (Bfield[2]/tesla)*1e9 << " nT" << std::endl;
}


