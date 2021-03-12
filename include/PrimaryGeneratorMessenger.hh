
// $Id: PrimaryGeneratorMessenger.hh,v 1.6 2002/12/16 16:37:26 maire Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
//

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;

class PrimaryGeneratorMessenger: public G4UImessenger
{
public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
  ~PrimaryGeneratorMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  PrimaryGeneratorAction*      Action;
  G4UIdirectory*               CRYDir;
  G4UIcmdWithAString*          FileCmd;
  G4UIcmdWithAString*          InputCmd;
  G4UIcmdWithoutParameter*     UpdateCmd;
  G4UIcmdWithADoubleAndUnit* cryDistanceZCmd;
  std::string* MessInput;
};

#endif
