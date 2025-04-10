

#include "SteppingActionMessenger.hh"

#include "SteppingAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIdirectory.hh"


SteppingActionMessenger::SteppingActionMessenger(SteppingAction* step)
  :G4UImessenger(),
   fSteppingAction(step) 
{
  fSteppingActionDir = new G4UIdirectory("/dataCollection/");
  fSteppingActionDir->SetGuidance("Select data collection type.");

  fcmd2 = new G4UIcmdWithAString("/dataCollection/setBackscatterFilename", this);
  fcmd2->SetParameterName("Enter file name.", true);
  fcmd2->SetDefaultValue("backscatter.csv");
}

SteppingActionMessenger::~SteppingActionMessenger()
{
  delete fSteppingActionDir;
  delete fcmd;
  delete fcmd2;
}

void SteppingActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fcmd2){
    fSteppingAction->SetBackscatterFilename(newValue);
  }
}
