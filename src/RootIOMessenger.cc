#include "RootIOMessenger.hh"
#include "RootIO.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
//#include <iostream.h>

RootIOMessenger::RootIOMessenger(RootIO * theIO)
:fIO(theIO)
{
    
    fRootDir = new G4UIdirectory("/root/");
    fRootDir->SetGuidance("Root IO control");
    
    theName = new G4UIcmdWithAString("/root/fileName",this);
    theName->SetGuidance("Set name of output root file.");
    theName->SetParameterName("fileName",true);
    theName->SetDefaultValue("hits.root");
    fIO->SetFileName( "hits.root" );
}

RootIOMessenger::~RootIOMessenger()
{
    delete fRootDir;
    delete theName;
}

void RootIOMessenger::SetNewValue( G4UIcommand * command, G4String newValues)
{
    if( command==theName ) { fIO->SetFileName(newValues);}
}

G4String RootIOMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    
    if( command==theName ){ cv = fIO->GetFileName(); }
    return cv;
}
