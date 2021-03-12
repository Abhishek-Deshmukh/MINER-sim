#include "PVolumeStore.hh"
#include <sstream>


#include "G4VPhysicalVolume.hh"

PVolumeStore::PVolumeStore(){}

PVolumeStore::~PVolumeStore(){}

void PVolumeStore::AddPVolume(const G4GeometryCell &cell){

    SetGeometryCell::iterator it =
    fSetGeometryCell.find(cell);
    if (it != fSetGeometryCell.end()) {
        G4cout << "PVolumeStore::AddPVolume: cell already stored"
        << G4endl;
        return;
    }

    fSetGeometryCell.insert(cell);


}

const G4VPhysicalVolume *PVolumeStore::
GetPVolume(const G4String &name) const {
    const G4VPhysicalVolume *pvol = 0;
    for (SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
         it != fSetGeometryCell.end(); ++it) {
        const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
        if (vol.GetName() == name) {
            pvol =  &vol;
        }
    }
    if (!pvol) {
        G4cout << "PVolumeStore::GetPVolume: no physical volume named: "
        << name << ", found" << G4endl;
    }
    return pvol;
}

G4String PVolumeStore::GetPNames() const {
    G4String NameString;
    for (SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
         it != fSetGeometryCell.end(); ++it) {
        const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
        std::ostringstream os;
        os << vol.GetName() << "_" << it->GetReplicaNumber()
        << "\n";
        G4String cellname = os.str();

        //    G4String cellname(vol.GetName());
        //    cellname += G4String("_");
        //    cellname += std::str(it->GetReplicaNumber());

        NameString += cellname;
    }
    return NameString;
}

G4int PVolumeStore::Size() {
    return fSetGeometryCell.size();
}
