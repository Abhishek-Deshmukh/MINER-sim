#ifndef MonitorGeometryConstruction_hh
#define MonitorGeometryConstruction_hh  MonitorGeometryConstruction_hh

#include "globals.hh"
#include <map>
#include <vector>
#include "G4GeometryCell.hh"
#include "PVolumeStore.hh"

#include "G4VUserParallelWorld.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class MonitorGeometryConstruction : public G4VUserParallelWorld
{
    public:
    MonitorGeometryConstruction(G4String worldName);
    virtual ~MonitorGeometryConstruction();
    
    G4VPhysicalVolume &GetWorldVolumeAddress() const;
    
    G4VPhysicalVolume* GetWorldVolume();
    
    void SetSensitive();
    
    
    virtual void ConstructSD();
    
    
    private:
    
    G4VPhysicalVolume* fGhostWorld;
    
};

#endif
