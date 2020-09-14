#include "globals.hh"
#include <sstream>

#include "ImportanceGeometryConstruction.hh"

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ImportanceGeometryConstruction::
ImportanceGeometryConstruction(G4String worldName)
:G4VUserParallelWorld(worldName),fLogicalVolumeVector()
{
    //  Construct();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ImportanceGeometryConstruction::~ImportanceGeometryConstruction()
{
    fLogicalVolumeVector.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ImportanceGeometryConstruction::SetSensitive(){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ImportanceGeometryConstruction::ConstructSD()
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


