

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"


RunActionMessenger::RunActionMessenger(RunAction* runAct)
  : G4UImessenger(),
    fRunAction(runAct) 
{
  fPrimDir = new G4UIdirectory("/dataCollection/");
  fPrimDir->SetGuidance("Set file names for histogram results files.");

  fcmd1 = new G4UIcmdWithAString("/dataCollection/setHistFileName",this);
  
  fcmd1->SetParameterName("Enter file name.",true);
  fcmd1->SetDefaultValue("energyHistogram.csv");
  fcmd1->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunActionMessenger::~RunActionMessenger()
{
  delete fPrimDir;
  delete fcmd1;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{
  if(command == fcmd1){
    fRunAction->SetHistFileName(newValue);
  }    	  

}
