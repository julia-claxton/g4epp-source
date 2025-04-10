

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"


RunActionMessenger::RunActionMessenger(RunAction* runAct)
  : G4UImessenger(),
    fRunAction(runAct) 
{
  fPrimDir = new G4UIdirectory("/dataCollection/");
  fPrimDir->SetGuidance("Set file names for energy deposition result file.");

  fcmd1 = new G4UIcmdWithAString("/dataCollection/setEnergyDepositionFileName", this);
  
  fcmd1->SetParameterName("Enter file name.",true);
  fcmd1->SetDefaultValue("energy_deposition.csv");
  fcmd1->AvailableForStates(G4State_PreInit, G4State_Idle);

  fcmd2 = new G4UIcmdWithAString("/dataCollection/setBackscatterFilename", this);
  fcmd2->SetParameterName("Enter file name.", true);
  fcmd2->SetDefaultValue("backscatter.csv");
}

RunActionMessenger::~RunActionMessenger()
{
  delete fPrimDir;
  delete fcmd1;
  delete fcmd2;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{
  if(command == fcmd1){
    fRunAction->SetEnergyDepositionFileName(newValue);
  }
  if(command == fcmd2){
    fRunAction->SetBackscatterFilename(newValue);
  }
}
