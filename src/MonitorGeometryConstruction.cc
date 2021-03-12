#include "globals.hh"
#include <sstream>

#include "MonitorGeometryConstruction.hh"
#include "MinerSD.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// For Primitive Scorers
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"

// for importance biasing
#include "G4IStore.hh"

MonitorGeometryConstruction::
MonitorGeometryConstruction(G4String worldName)
:G4VUserParallelWorld(worldName)
{
    //  Construct();
}

MonitorGeometryConstruction::~MonitorGeometryConstruction()
{
}

void MonitorGeometryConstruction::SetSensitive(){

}
void MonitorGeometryConstruction::ConstructSD()
{

    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    MinerSD* insideIB_Det = new MinerSD("/MINERsim/insideIB_Det","MS_insideIB_hits");
    SDman->AddNewDetector(insideIB_Det);
    SetSensitiveDetector("insideIB_log",insideIB_Det, true);

    MinerSD* atMovable_Det = new MinerSD("/MINERsim/atMovable_Det","MS_atMovable_hits");
    SDman->AddNewDetector(atMovable_Det);
    SetSensitiveDetector("atMovable_log",atMovable_Det, true);

    MinerSD* atNeutronDet_Det = new MinerSD("/MINERsim/atNeutronDet_Det","MS_atNeutronDet_hits");
    SDman->AddNewDetector(atNeutronDet_Det);
    SetSensitiveDetector("atNeutronDet_log",atNeutronDet_Det, true);

    MinerSD* atLead_Det = new MinerSD("/MINERsim/atLead_Det","MS_atLead_hits");
    SDman->AddNewDetector(atLead_Det);
    SetSensitiveDetector("atLead_log",atLead_Det, true);

    MinerSD* atStartTC_Det = new MinerSD("/MINERsim/atStartTC_Det","MS_atStartTC_hits");
    SDman->AddNewDetector(atStartTC_Det);
    SetSensitiveDetector("atStartTC_log",atStartTC_Det, true);
}
