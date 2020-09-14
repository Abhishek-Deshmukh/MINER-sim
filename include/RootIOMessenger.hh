#ifndef RootIOMessenger_h
#define RootIOMessenger_h 1

class RootIO;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class RootIOMessenger: public G4UImessenger
{
public:
    RootIOMessenger(RootIO * theIO);
    ~RootIOMessenger();
    
public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
    
private:
    RootIO * fIO;
    
private: //commands
    G4UIdirectory*              fRootDir;
    G4UIcmdWithAString *        theName;
    
};

#endif
