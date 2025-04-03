

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class PrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;

class PrimaryGeneratorMessenger : public G4UImessenger
{
public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction* );
  virtual ~PrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);


private:
  PrimaryGeneratorAction* fPrimaryGenerator;
  G4UIdirectory*             fPrimDir;
  G4UIcmdWithAnInteger*      fcmd;
  G4UIcmdWithAnInteger*      fcmd2;
  G4UIcmdWithAnInteger*      fcmd3;
  G4UIcmdWithADouble*        fDcmd;
  G4UIcmdWithADouble*        fDcmd2;
  G4UIcmdWithADouble*        fDcmd3;

};

#endif
