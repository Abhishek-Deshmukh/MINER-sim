#ifndef ImportanceGeometryConstruction_hh
#define ImportanceGeometryConstruction_hh  ImportanceGeometryConstruction_hh

#include "globals.hh"
#include <map>
#include <vector>
#include "G4GeometryCell.hh"
#include "PVolumeStore.hh"

#include "G4VUserParallelWorld.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VIStore;
class G4VWeightWindowStore;



class ImportanceGeometryConstruction : public G4VUserParallelWorld
{
    public:
    ImportanceGeometryConstruction(G4String worldName);
    virtual ~ImportanceGeometryConstruction();
    
    const G4VPhysicalVolume &GetPhysicalVolumeByName(const G4String& name) const;
    G4VPhysicalVolume &GetWorldVolumeAddress() const;
    G4String ListPhysNamesAsG4String();
    G4String GetCellName(G4int i);
    //G4GeometryCell GetGeometryCell(G4int i);
    
    G4VPhysicalVolume* GetWorldVolume();
    
    void SetSensitive();
    
    virtual void Construct();
    virtual void ConstructSD();
    
    G4VIStore* CreateImportanceStore();
    // create an importance store, caller is responsible for deleting it
    
    
    private:
    PVolumeStore fPVolumeStore;
    
    
    
    //  std::vector< G4VPhysicalVolume * > fPhysicalVolumeVector;
    std::vector< G4LogicalVolume * > fLogicalVolumeVector;
    
    //  G4VPhysicalVolume *fWorldVolume;
    
   // G4VPhysicalVolume* fGhostWorld;
    
};

#endif
