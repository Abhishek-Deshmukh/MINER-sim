//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//


#ifndef GeometryConstruction1_h
#define GeometryConstruction1_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4GDMLParser.hh"


class G4VPhysicalVolume;
class G4Material;
class G4VIStore;
class G4VWeightWindowStore;
class DetectorMessenger;

class GeometryConstruction : public G4VUserDetectorConstruction
{
public:
    GeometryConstruction ();
    ~GeometryConstruction ();



    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    G4VPhysicalVolume* ConstructDetector();
    G4LogicalVolume* Detector();
    G4LogicalVolume* BaseConstruct();
    G4LogicalVolume* PentagonConstruct();



    G4VIStore* CreateImportanceStore();
    // create an importance store, caller is responsible for deleting it

    G4VWeightWindowStore *CreateWeightWindowStore();
    // create an weight window  store, caller is responsible for
    // deleting it


    void SetSensitive();
    void ListOfMaterials();

    void SetReadFile( const G4String& File );

    inline G4VIStore* GetGeomStore(){return aIstore;};


    G4VPhysicalVolume* GetWorldVolume();
    G4LogicalVolume* ConstructBEGe(std::string name);


private:
    //G4VPhysicalVolume* fUniverse_phys;
    G4Material* fAir ;
    G4Material* fAluminum ;
    G4Material* fPb;
    G4Material* fXenon;
    G4Material* fCu;
    G4Material* ssteel;
    G4Material* LN2;
    G4Material* fGe;
    G4Material* fCd;

    G4GDMLParser fParser;
    //std::vector< G4LogicalVolume * > fLogicalVolumeVector;
    //std::vector< G4VPhysicalVolume * > fPhysicalVolumeVector;

    //const MINERMaterials * mats = MINERMaterials::GetInstance();
    G4double in = 2.54*CLHEP::cm;
    G4bool fCheckOverlaps = true;
    G4RotationMatrix *zeroRot = new G4RotationMatrix;
    G4ThreeVector zeroPos;

    G4String fReadFile;
    G4int fWritingChoice;

    DetectorMessenger* fDetectorMessenger;


    G4VisAttributes* aVisAttWater = new G4VisAttributes(G4Colour(0.0,0,1.0));
    G4VisAttributes* aVisAttGraphite = new G4VisAttributes(G4Colour(.2,0.2,.2));
    G4VisAttributes* aVisAttLead = new G4VisAttributes(G4Colour(.1,0.1,.1));
    G4VisAttributes* aVisAttPoly = new G4VisAttributes(G4Colour(31./255.,191./255.,127./255.));
    G4VisAttributes* aVisAttAlTubes = new G4VisAttributes(G4Colour(.5,.5,.5));
    G4VisAttributes* aVisAttStainless = new G4VisAttributes(G4Colour(.5,.5,.5));
    G4VisAttributes* aVisAttHDConcrete = new G4VisAttributes(G4Colour(0.3,0.3,.3));
    G4VisAttributes* aVisAttBioShield = new G4VisAttributes(G4Colour(0.55,0.55,.55));
    G4VisAttributes* aVisAttSSLiner = new G4VisAttributes(G4Colour(0.4,0.4,.4));
    G4VisAttributes* aVisAttLiqN2 = new G4VisAttributes(G4Colour(0.0,0.0,.8));
    G4VisAttributes* aVisAttDet = new G4VisAttributes(G4Colour(145./255.,151./255.,171./255.));
    G4VisAttributes* aVisAttCu = new G4VisAttributes(G4Colour(252./255.,187./255.,34./255.));
    G4VisAttributes* aVisAttScint = new G4VisAttributes(G4Colour(245./255.,241./255.,12./255.));
    //G4VisAttributes* aVisAttWhite = new G4VisAttributes(G4Colour::White());

    G4VIStore *aIstore;


};

#endif
