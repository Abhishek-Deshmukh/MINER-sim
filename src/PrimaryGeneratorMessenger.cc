//
// $Id: PrimaryGeneratorMessenger.cc,v 1.8 2002/12/16 16:37:27 maire Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
//

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
  :Action(Gun)
{
  CRYDir = new G4UIdirectory("/CRY/");
  CRYDir->SetGuidance("CRY initialization");

  FileCmd = new G4UIcmdWithAString("/CRY/file",this);
  FileCmd->SetGuidance("This reads the CRY definition from a file");
  FileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  InputCmd = new G4UIcmdWithAString("/CRY/input",this);
  InputCmd->SetGuidance("CRY input lines");
  InputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/CRY/update",this);
  UpdateCmd->SetGuidance("Update CRY definition.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed the CRY definition.");
  UpdateCmd->AvailableForStates(G4State_Idle);

  cryDistanceZCmd = new G4UIcmdWithADoubleAndUnit("/CRY/distanceZ", this);
  cryDistanceZCmd->SetGuidance("Set Z distance of the plane from which the cosmic rays are generated.");
  cryDistanceZCmd->SetParameterName("cryDistanceZCmd", false);
  cryDistanceZCmd->SetDefaultUnit("m");

  MessInput = new std::string;

}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete CRYDir;
  delete InputCmd;
  delete UpdateCmd;
  delete FileCmd;
  delete cryDistanceZCmd;
}

void PrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{
  // if( command == InputCmd )
  //  {
  //    Action->InputCRY();
  //    (*MessInput).append(newValue);
  //    (*MessInput).append(" ");
  //  }

  if( command == UpdateCmd )
   {
     Action->UpdateCRY(MessInput);
     *MessInput = "";
   }

  if( command == FileCmd )
   { Action->UpdateCRYFile(newValue); }

  if (command == cryDistanceZCmd) {
    Action->SetSourceZPosition(cryDistanceZCmd->GetNewDoubleValue(newValue));
  }

}
