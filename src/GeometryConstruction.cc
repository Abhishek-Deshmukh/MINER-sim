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
// $Id: GeometryConstruction.cc 67272 2013-02-13 12:05:47Z ihrivnac $
//
/// \file eventgenerator/exgps/src/GeometryConstruction.cc
/// \brief Implementation of the GeometryConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "GeometryConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "MinerSD.hh"
//#include "MINERMaterials.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"

#include "G4GDMLParser.hh"


//#include "DetectorMessenger.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GeometryConstruction::GeometryConstruction()
: G4VUserDetectorConstruction(),
fAir(0), fAluminum(0), fPb(0), fXenon(0), fCu(0), ssteel(0), LN2(0), fGe(0), fCd(0)
{
    
    fReadFile ="MINERSetupAugust30.gdml";
//    fReadFile ="MinerDet.gdml";
    fWritingChoice=1;
    
    //fDetectorMessenger = new DetectorMessenger( this );
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GeometryConstruction::~GeometryConstruction()
{
    
    //if(fDetectorMessenger) delete fDetectorMessenger;
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GeometryConstruction::ConstructSDandField()
{
    
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    // Sensitive detectors
    //MinerSD* aCraigSD = new MinerSD("veloSD", "MS_Craig_hits");
    
    MinerSD* HybridSD1 = new MinerSD("HybridSD", "Hybrid_hits");
    MinerSD* HVSD1 = new MinerSD("HVSD", "HV_hits");
    MinerSD* CsISD1 = new MinerSD("CsISD", "CsI_hits");

    //SDman->AddNewDetector(aCraigSD);
    
    SDman->AddNewDetector(HybridSD1);
    SDman->AddNewDetector(HVSD1);
    SDman->AddNewDetector(CsISD1);
    
    
    //SetSensitiveDetector("fLogicDet", aCraigSD, true);
    
   // MinerSD* aCraigSD2 = new MinerSD("/MINERsim/Craig2", "MS_Craig_hits2");
    //SDman->AddNewDetector(aCraigSD2);
   // SetSensitiveDetector("Craig_log_outsideCastle", aCraigSD2, true);
    
    
    const G4GDMLAuxMapType* auxmap = fParser.GetAuxMap();
    G4cout << "Found " << auxmap->size()
    << " volume(s) with auxiliary information."
    << G4endl << G4endl;
    for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
        iter!=auxmap->end(); iter++)
    {
        G4cout << "Volume " << ((*iter).first)->GetName()
        << " has the following list of auxiliary information: "
        << G4endl << G4endl;
        for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
             vit!=(*iter).second.end(); vit++)
        {
            G4cout << "--> Type: " << (*vit).type
            << " Value: " << (*vit).value << G4endl;
        }
    }
    G4cout << G4endl;
    
    // The same as above, but now we are looking for
    // sensitive detectors setting them for the volumes
    
    for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
        iter!=auxmap->end(); iter++)
    {
        G4cout << "Volume " << ((*iter).first)->GetName()
        << " has the following list of auxiliary information: "
        << G4endl << G4endl;
        for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
             vit!=(*iter).second.end();vit++)
        {
            if ((*vit).type=="SensDet")
            {
                G4cout << "Attaching sensitive detector " << (*vit).value
                << " to volume " << ((*iter).first)->GetName()
                <<  G4endl << G4endl;
                
                G4VSensitiveDetector* mydet =
                SDman->FindSensitiveDetector((*vit).value);
                if(mydet)
                {
                    G4LogicalVolume* myvol = (*iter).first;
                    myvol->SetSensitiveDetector(mydet);
                }
                else
                {
                    G4cout << (*vit).value << " detector not found" << G4endl;
                }
            }
        }
    }
    
    
}

G4VPhysicalVolume* GeometryConstruction::Construct()
{
    
    G4VPhysicalVolume* fUniverse_phys;
    
   
        fParser.Read(fReadFile);
      
        G4cout << *(G4Material::GetMaterialTable() ) << G4endl;
        
        // Giving World Physical Volume from GDML Parser
        //
        fUniverse_phys = fParser.GetWorldVolume();
    
    
    G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    fUniverse_phys->GetLogicalVolume()->SetVisAttributes(BoxVisAtt);
    
    return fUniverse_phys;
}

void GeometryConstruction::ListOfMaterials()
{
   /* G4double a;  // atomic mass
    G4double z;  // atomic number
    G4double density,temperature,pressure;
    G4double fractionmass;
    G4String name, symbol;
    G4int ncomponents;
    
    // Elements needed for the materials
    
    a = 14.01*g/mole;
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);
    
    a = 16.00*g/mole;
    G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);
    
    a = 26.98*g/mole;
    G4Element* elAl = new G4Element(name="Aluminum", symbol="Al", z=13., a);
    
    a=55.85*g/mole;
    G4Element* Fe = new G4Element(name="Iron" ,symbol="Fe" , z= 26., a);
    
    a=12.011*g/mole;
    G4Element* C  = new G4Element( name="Carbon", symbol="C", z= 6. , a);
    
    a=58.9332*g/mole;
    G4Element* Co = new G4Element( name="Cobalt", symbol="Co", z=27. , a);
    
    // Print the Element information
    //
    G4cout << *(G4Element::GetElementTable()) << G4endl;
    
    //SSteel
    
    density = 7.7 *g/cm3;
    ssteel = new G4Material(name="ssteel", density, ncomponents=3);
    ssteel->AddElement(C, fractionmass = 0.04);
    ssteel->AddElement(Fe, fractionmass =0.88);
    ssteel->AddElement(Co, fractionmass =0.08);
    
    
    LN2 = new G4Material("LN2", z=7., a= 14.00*g/mole, density=0.807*g/cm3, kStateLiquid, temperature= 77.*kelvin, pressure= 1.0*atmosphere);
    
    
    // Air
    //
    density = 1.29*mg/cm3;
    fAir = new G4Material(name="Air", density, ncomponents=2);
    fAir->AddElement(elN, fractionmass=0.7);
    fAir->AddElement(elO, fractionmass=0.3);
    
    // Aluminum
    //
    density = 2.70*g/cm3;
    fAluminum = new G4Material(name="Aluminum", density, ncomponents=1);
    fAluminum->AddElement(elAl, fractionmass=1.0);
    
    // Lead
    //
    fPb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
    
    fCu = new G4Material("Copper", z=29., a=63.55*g/mole, density = 8.96 *g/cm3);
    
    fGe = new G4Material("Germanium", z=32., a=72.63*g/mole, density = 5.323 *g/cm3);
    
    fCd = new G4Material("Cadmium", z=48., a=112.414*g/mole, density = 8.65 *g/cm3);
    
    // Xenon gas
    //
    fXenon = new G4Material("XenonGas", z=54., a=131.29*g/mole,
                            density= 5.458*mg/cm3, kStateGas,
                            temperature= 293.15*kelvin, pressure= 1*atmosphere);
    
    // Prints the material information
    //
    G4cout << *(G4Material::GetMaterialTable() ) << G4endl;*/
}



G4VPhysicalVolume* GeometryConstruction::ConstructDetector()
{
    /*const G4double fExpHall_x = 5.*m;
    // Arbitary values that should enclose any reasonable geometry
    //
    const G4double expHall_y = fExpHall_x;
    const G4double expHall_z = fExpHall_x;
    
    G4RotationMatrix * rotX = new G4RotationMatrix();
    rotX->rotateX(270*deg);
    
    G4RotationMatrix * rotY = new G4RotationMatrix();
    rotY->rotateY(90*deg);
    
    
    //G4LogicalVolume *detectorLV = ConstructBEGe("_insideCastle");
    ///WORLD
    
    
    G4Box * experimentalHallBox =
    new G4Box("ExpHallBox", fExpHall_x, expHall_y, expHall_z);
    G4LogicalVolume * experimentalHallLV =
    new G4LogicalVolume(experimentalHallBox, fAir, "ExpHallLV");
    G4PVPlacement * experimentalHallPhys =
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), experimentalHallLV,
                      "ExpHallPhys", 0, false, 0);
    
    
    
    G4LogicalVolume* detector = ConstructBEGe("_insideCastle");
    new G4PVPlacement(rotX, G4ThreeVector(0.0,100.0,0.0), detector, "DetConst", experimentalHallLV, false, 0);
    
    
    
    G4LogicalVolume* baseShielding = BaseConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,-202.25,0.0), baseShielding, "baseConst", experimentalHallLV, false, 0);
    G4LogicalVolume* baseShielding2 = BaseConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,-304.0,0.0), baseShielding2, "baseConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,-100.5*mm,0.0), pentagonShielding, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding1 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1.25*mm,0.0), pentagonShielding1, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding2 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,103.0*mm,0.0), pentagonShielding2, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding3 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,204.75*mm,0.0), pentagonShielding3, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding4 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,306.50*mm,0.0), pentagonShielding4, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding5 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,408.25*mm,0.0), pentagonShielding5, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding6 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,510.5*mm,0.0), pentagonShielding6, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding7 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,612.25*mm,0.0), pentagonShielding7, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding8 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,714.0*mm,0.0), pentagonShielding8, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding9 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,815.75*mm,0.0), pentagonShielding9, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding10 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,917.5*mm,0.0), pentagonShielding10, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding11 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1019.25*mm,0.0), pentagonShielding11, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding12 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1121.0*mm,0.0), pentagonShielding12, "pentConst", experimentalHallLV, false, 0);
    G4LogicalVolume* pentagonShielding13 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1222.5*mm,0.0), pentagonShielding13, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding14 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1324.5*mm,0.0), pentagonShielding14, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding15 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1426.25*mm,0.0), pentagonShielding15, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding16 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1528.0*mm,0.0), pentagonShielding16, "pentConst", experimentalHallLV, false, 0);
    
    G4LogicalVolume* pentagonShielding17 = PentagonConstruct();
    new G4PVPlacement(rotY, G4ThreeVector(0.0,1629.75*mm,0.0), pentagonShielding17, "pentConst", experimentalHallLV, false, 0);
    
    return experimentalHallPhys;*/
}


G4LogicalVolume*GeometryConstruction::ConstructBEGe(std::string name)
{
    
    /*G4Box * detectorBox =
    new G4Box("detectorBox", 10000.0, 10000.0, 10000.0);
    G4LogicalVolume * detectorLV =
    new G4LogicalVolume(detectorBox, fAir, "detLV");
    // G4PVPlacement * detectorPhys =
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), detectorLV,
                      "detPhys", 0, false, 0);
    
    
    
    G4double detRad = 50.8*mm/2.;
    G4double detHalfZ = (4.805/2.)*cm;
    
    G4VSolid* target1 = new G4Tubs("Target1"+name,0.,detRad,detHalfZ,0,360*deg);
    G4VSolid* target2 = new G4Tubs("Target2"+name,0.,4.75*mm,15.0*mm,0,360*deg);
    G4ThreeVector targetOffset(0,0,-1*(detHalfZ - 15.0*mm));
    G4SubtractionSolid *target = new G4SubtractionSolid("Target"+name,target1,target2,0,targetOffset);
    G4LogicalVolume *fLogicDet = new G4LogicalVolume(target, fGe,"Craig_log"+name);
    
    
    G4double cuCasingThick = 0.76*mm;
    //G4double cuCasingThick = 2.25*mm;
    G4double cuRingThick = 2.7*mm - 0.76*mm;
    //G4double cuRingThick = 3.5*mm;
    G4double cuRingHeight = 8.6*mm;
    
    
    G4VSolid* cuCasing1 = new G4Tubs("cuCasing1"+name,0.,detRad+cuCasingThick,35.5*mm/2. + detHalfZ + cuCasingThick/2.,0,360*deg);
    G4VSolid* cuCasing2 = new G4Tubs("cuCasing2"+name,0.,detRad+.001*mm,35.5*mm/2. + detHalfZ + cuCasingThick/2.,0,360*deg);
    G4ThreeVector casingOffset(0,0,3.0*mm);
    G4SubtractionSolid *cuCasing = new G4SubtractionSolid("cuCasting"+name,cuCasing1,cuCasing2,0,casingOffset);
    
    
    G4VSolid* cuRing = new G4Tubs("cuRing"+name,detRad+cuCasingThick,detRad+cuCasingThick+cuRingThick,cuRingHeight/2.,0,360*deg);
    G4VSolid* cuRingT = new G4Tubs("cuRingT"+name,detRad+cuCasingThick,detRad+cuCasingThick+cuRingThick,cuRingHeight/6.,0,360*deg);
    G4ThreeVector ringOffsetT(0,0,35.5*mm/2. + detHalfZ + cuCasingThick/2.-cuRingHeight/6.);
    G4ThreeVector ringOffset1(0,0,(35.5*mm/2. + detHalfZ + cuCasingThick/2.)/3.);
    G4ThreeVector ringOffset2(0,0,-1*(35.5*mm/2. + detHalfZ + cuCasingThick/2.)/3.);
    
    G4UnionSolid *cuDetHousing1 = new G4UnionSolid("cuDetHousing1"+name,cuCasing,cuRingT,0,ringOffsetT);
    G4UnionSolid *cuDetHousing2 = new G4UnionSolid("cuDetHousing2"+name,cuDetHousing1,cuRing,0,ringOffset1);
    G4UnionSolid *cuDetHousing = new G4UnionSolid("cuDetHousing"+name,cuDetHousing2,cuRing,0,ringOffset2);
    
    G4LogicalVolume *fLogicCasing = new G4LogicalVolume(cuDetHousing2, fCu,"DetHousing_log"+name);
    G4ThreeVector cuHousingPos(0,0,-460 + (-1*((35.5*mm/2. + detHalfZ + cuCasingThick/2.) - detHalfZ)));
    //G4PVPlacement * fLogicCasingPos =
    new G4PVPlacement(0, cuHousingPos, fLogicCasing, "CasingPos", detectorLV, false, 0);
    
    
    G4double tubeAthick = (1.59)*mm;
    G4double tubeAlen = 111.25*mm;
    G4double tubeAdia = 76.2*mm;
    
    G4double tubeBthick = (3.18)*mm;
    G4double tubeBlen = 127*mm;
    G4double tubeBdia = 79.5*mm;
    
    G4double tubeCthick = (3.18)*mm;
    G4double tubeClen = 22.35*mm;
    G4double tubeCdia = 28.7*mm;
    
    G4double tubeDthick = (3.18)*mm;
    G4double tubeDlen = 349.25*mm;
    G4double tubeDdia = 228.6*mm;
    
    G4Tubs *outTubeA = new G4Tubs("outTubeA"+name, 0., tubeAdia/2., tubeAlen/2.,0,360*deg);
    G4Tubs *inTubeA = new G4Tubs("inTubeA"+name, 0., tubeAdia/2. - tubeAthick, tubeAlen/2. - tubeAthick/2.,0,360*deg);
    // right now the above defines the aluminum tube surrounding the detector, and has a 1/16" aluminum window, probably too thick?
    G4ThreeVector inTubeAOffset(0,0,-tubeAthick/2.);
    G4SubtractionSolid *detTubeA = new G4SubtractionSolid("detTubeA"+name,outTubeA,inTubeA,0,inTubeAOffset);
    
    G4Tubs *outTubeB = new G4Tubs("outTubeB"+name, 0., tubeBdia/2., tubeBlen/2.,0,360*deg);
    G4Tubs *inTubeB = new G4Tubs("inTubeB"+name, 0., tubeBdia/2. - tubeBthick, tubeBlen/2. - tubeBthick/2.,0,360*deg);
    G4ThreeVector inTubeBOffset(0,0,+tubeBthick/2.);
    G4SubtractionSolid *detTubeB = new G4SubtractionSolid("detTubeB"+name,outTubeB,inTubeB,0,inTubeBOffset);
    
    G4Tubs *detTubeC = new G4Tubs("detTubeC"+name, tubeCdia/2. - tubeCthick, tubeCdia/2., tubeClen/2.,0,360*deg);
    //G4Tubs *inTubeC = new G4Tubs("inTubeC"+name, 0., tubeCdia/2. - tubeCthick, tubeClen/2.,0,360*deg);
    
    G4Tubs *outTubeD = new G4Tubs("outTubeD"+name, 0., tubeDdia/2., tubeDlen/2.,0,360*deg);
    G4Tubs *inTubeD = new G4Tubs("inTubeD"+name, 0., tubeDdia/2. - tubeDthick, tubeDlen/2. - tubeDthick,0,360*deg);
    G4SubtractionSolid *detTubeD = new G4SubtractionSolid("detTubeD"+name,outTubeD,inTubeD);
    
    
    
    G4ThreeVector tubeAPos(0.,0.,               -350 + ( -1*(tubeAlen/2. - tubeAthick - detHalfZ) + 5*mm)); //detector 5mm from front aluminum
    G4ThreeVector tubeBPos(0,0,                -350 + ((-1*(tubeAlen/2. - tubeAthick - detHalfZ) + 5*mm) - tubeAlen/2. - tubeBlen/2.));
    G4ThreeVector tubeCPos(0,0 - (tubeBdia/3.), ((-1*(tubeAlen/2. - tubeAthick - detHalfZ) + 5*mm) - tubeAlen/2. - tubeBlen - tubeClen/2.));
    G4ThreeVector tubeDPos(0,0,                 ((-1*(tubeAlen/2. - tubeAthick - detHalfZ) + 5*mm) - tubeAlen/2. - tubeBlen - tubeClen - tubeDlen/2.));
    
    G4ThreeVector tubeA_N2_Pos(0, 0 , -25.805);//25.805
    G4ThreeVector tubeB_N2_Pos (0, 0, ((-1*(tubeAlen/2. - tubeAthick - detHalfZ) + 5*mm) - tubeAlen/2. - tubeBlen/2.) +tubeBthick/2.) ;
    
    G4Tubs *liqN2_TubeD = new G4Tubs("liqN2_TubeD"+name, 0., tubeDdia/2. - tubeDthick, tubeDlen/2. - tubeDthick,0,360*deg);
    
    // canberra lead shield
    G4Tubs *InnerCanberra = new G4Tubs("InnerCanberra"+name, 0., (279.*mm)/2.,203.*mm,0.,360*deg);
    G4Tubs *CopperCanberraI = new G4Tubs("CopperCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm,203.*mm + 1.6*mm,0.,360*deg);
    G4Tubs *TinCanberraI = new G4Tubs("TinCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm + 1.0*mm,203.*mm + 1.6*mm + 1.0*mm,0.,360*deg);
    G4Tubs *LeadCanberraI = new G4Tubs("LeadCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm + 1.0*mm + 101.6*mm,203.*mm + 1.6*mm + 1.0*mm + 101.6*mm,0.,360*deg);
    G4Tubs *SteelCanberraI = new G4Tubs("SteelCanberraI"+name, 0., (279.*mm)/2. + 1.6*mm + 1.0*mm + 101.6*mm + 9.5*mm,203.*mm + 1.6*mm + 1.0*mm + 101.6*mm + 9.5*mm,0.,360*deg);
    
    G4SubtractionSolid *CopperCanberraI2 = new G4SubtractionSolid("CopperCanberraI2"+name,CopperCanberraI,InnerCanberra);
    G4SubtractionSolid *TinCanberraI2 = new G4SubtractionSolid("TinCanberraI2"+name,TinCanberraI,CopperCanberraI);
    G4SubtractionSolid *LeadCanberraI2 = new G4SubtractionSolid("LeadCanberraI2"+name,LeadCanberraI,TinCanberraI);
    G4SubtractionSolid *SteelCanberraI2 = new G4SubtractionSolid("SteelCanberraI2"+name,SteelCanberraI,LeadCanberraI);
    
    G4Tubs *bottomHoleCanberra = new G4Tubs("bottomHoleCanberra"+name,0.,(3.25)*mm/2.,(1.6*mm + 1.0*mm + 101.6*mm + 9.5*mm)/2.,0.,360*deg);
    G4ThreeVector bHC_Pos(0,0,-1*(1.6*mm + 1.0*mm + 101.6*mm + 9.5*mm)/2. - 203.*mm);
    
    G4SubtractionSolid *CopperCanberra = new G4SubtractionSolid("CopperCanberra"+name,CopperCanberraI2,bottomHoleCanberra,0,bHC_Pos);
    G4SubtractionSolid *TinCanberra = new G4SubtractionSolid("TinCanberra"+name,TinCanberraI2,bottomHoleCanberra,0,bHC_Pos);
    G4SubtractionSolid *LeadCanberra = new G4SubtractionSolid("LeadCanberra"+name,LeadCanberraI2,bottomHoleCanberra,0,bHC_Pos);
    G4SubtractionSolid *SteelCanberra = new G4SubtractionSolid("SteelCanberra"+name,SteelCanberraI2,bottomHoleCanberra,0,bHC_Pos);
    
    G4ThreeVector ShieldCanberra_Pos(0,0,17 + (645*mm/2. + tubeDPos.z() + tubeDlen/2.));
    
    
    G4LogicalVolume *fLogicDetTubeA = new G4LogicalVolume(detTubeA, fAluminum,"detTubeA_log"+name);
    G4LogicalVolume *fLogicDetTubeB = new G4LogicalVolume(detTubeB, fAluminum,"detTubeB_log"+name);
    G4LogicalVolume *fLogicDetTubeC = new G4LogicalVolume(detTubeC, ssteel,"detTubeC_log"+name);
    G4LogicalVolume *fLogicDetTubeD = new G4LogicalVolume(detTubeD, ssteel,"detTubeD_log"+name);
    G4LogicalVolume *fLogicDetLiN2D = new G4LogicalVolume(liqN2_TubeD,LN2,"liqN2_TubeD_log"+name);
    //G4LogicalVolume *fLogicDetLead = new G4LogicalVolume(detLeadShield,fAluminum,"detLeadShield_log"+name);
    G4LogicalVolume *fLogicDetCSCu = new G4LogicalVolume(CopperCanberra,fCu,"CopperCanberra_log"+name);
    //G4LogicalVolume *fLogicDetCSTin = new G4LogicalVolume(TinCanberra,fAluminum,"TinCanberra_log"+name);
    G4LogicalVolume *fLogicDetCSTin = new G4LogicalVolume(TinCanberra,fCd,"TinCanberra_log"+name);
    G4LogicalVolume *fLogicDetCSPb = new G4LogicalVolume(LeadCanberra,fPb,"LeadCanberra_log"+name);
    G4LogicalVolume *fLogicDetCSSteel = new G4LogicalVolume(SteelCanberra,ssteel,"SteelCanberra_log"+name);
    
    
    new G4PVPlacement(0, tubeAPos,fLogicDetTubeA, "TubeAPos", detectorLV, false, 0);
    
    new G4PVPlacement(0, tubeBPos,fLogicDetTubeB, "TubeBPos", detectorLV, false, 0);
    //G4PVPlacement * fLogicDetTubeCPos =
    new G4PVPlacement(0,tubeCPos, fLogicDetTubeC, "TubeCPos",detectorLV, false, 0);
    //G4PVPlacement * fLogicDetTubeDPos =
    new G4PVPlacement(0, tubeDPos, fLogicDetTubeD,"TubeDPos", detectorLV, false, 0);
    //G4PVPlacement * fLogicDetTubeLiN2DPos =
    new G4PVPlacement(0,tubeDPos, fLogicDetLiN2D,"TubeLiN2DPos", detectorLV, false, 0);
    //G4PVPlacement * fLogicPet =
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -460.0), fLogicDet , "fLogicDetPos", detectorLV, false, 0);
    
    //fullDet->AddPlacedVolume(fLogicDetLead,detLead_Pos,zeroRot);
    //G4PVPlacement * fLogicDetCSCuPos =
    new G4PVPlacement(0, ShieldCanberra_Pos, fLogicDetCSCu, "DetCSCuPos", detectorLV, false, 0);
    //G4PVPlacement * fLogicDetCSTinPos =
    new G4PVPlacement(0, ShieldCanberra_Pos, fLogicDetCSTin, "DetCSTinPos", detectorLV, false, 0);
    //G4PVPlacement * fLogicDetCSPbPos =
    new G4PVPlacement(0, ShieldCanberra_Pos, fLogicDetCSPb, "DetCSPbPos", detectorLV, false, 0);
    //G4PVPlacement * fLogicDetCSSteelPos =
    new G4PVPlacement(0, ShieldCanberra_Pos, fLogicDetCSSteel, "DetCSSteelPos", detectorLV, false, 0);
    
    
    return detectorLV;
    
    */
    
}
G4LogicalVolume*GeometryConstruction::BaseConstruct()
{
    /*const G4double br_x = 10.16*cm;
    const G4double br_y = 5.08*cm;
    const G4double br_z = 20.32*cm;
    const G4double bx_l = 30.48*cm;
    const G4double rct_l = 50.80*cm;
    
    
    
    G4Box * detectorBox =
    new G4Box("detectorBox", 10000.0, 10000.0, 10000.0);
    G4LogicalVolume * detectorLV =
    new G4LogicalVolume(detectorBox, fAir, "detLV");
    // G4PVPlacement * detectorPhys =
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), detectorLV,
                      "detPhys", 0, false, 0);
    ///Center Square Brick
    G4Box * S_brck =
    new G4Box("Square", bx_l, br_y, bx_l);
    G4LogicalVolume * Square_br =
    new G4LogicalVolume(S_brck, fPb, "Square1");
    new G4PVPlacement(0, G4ThreeVector(-0.5, -200.0, 654.5), Square_br, "Sqr", detectorLV, false, 0);
    
    /////Six Bricks above Square
    G4VSolid * Brcks =
    new G4Box("Brcks", br_x, br_y, br_z);
    G4LogicalVolume * Brcks1=
    new G4LogicalVolume(Brcks, fPb, "BrcksLV");
    new G4PVPlacement(0, G4ThreeVector(203.5,-200.0*mm,146.0), Brcks1, "Brck1", detectorLV, false, 0);
    G4LogicalVolume * Brcks2 =
    new G4LogicalVolume(Brcks, fPb, "BrcksLV2");
    new G4PVPlacement(0, G4ThreeVector(0.0*mm, -200.0*mm, 146.0), Brcks2, "Brck2", detectorLV, false, 0);
    G4LogicalVolume * Brcks3=
    new G4LogicalVolume(Brcks, fPb, "BrcksLV3");
    new G4PVPlacement(0, G4ThreeVector(-203.5,-200.0*mm,146.0), Brcks3, "Brck3", detectorLV, false, 0);
    G4LogicalVolume * Brcks4 =
    new G4LogicalVolume(Brcks, fPb, "BrcksLV4");
    new G4PVPlacement(0, G4ThreeVector(-407.0*mm, -200.0*mm, 146.0), Brcks4, "Brck4", detectorLV, false, 0);
    G4LogicalVolume * Brcks5=
    new G4LogicalVolume(Brcks, fPb, "BrcksLV5");
    new G4PVPlacement(0, G4ThreeVector(-610.5,-200.0*mm, 146.0), Brcks5, "Brck5", detectorLV, false, 0);
    G4LogicalVolume * Brcks6 =
    new G4LogicalVolume(Brcks, fPb, "BrcksLV6");
    new G4PVPlacement(0, G4ThreeVector(-814.0*mm, -200.0*mm, 146.0), Brcks6, "Brck6", detectorLV, false, 0);
    ////Left of Square Bricks, Horizontal bricks above square
    G4Box * BrcksSW =
    new G4Box("Brcks2", br_z, br_y, br_x);
    
    G4LogicalVolume * Brcks7=
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV7");
    new G4PVPlacement(0, G4ThreeVector(-508.5,-200.0,451), Brcks7, "Brck7", detectorLV, false, 0);
    G4LogicalVolume * Brcks8 =
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV8");
    new G4PVPlacement(0, G4ThreeVector(-508.5*mm, -200.0, 654.5), Brcks8, "Brck8", detectorLV, false, 0);
    G4LogicalVolume * Brcks9=
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV9");
    new G4PVPlacement(0, G4ThreeVector(-508.5, -200.0, 858.0), Brcks9, "Brck9", detectorLV, false, 0);
    G4LogicalVolume * Brcks10 =
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV10");
    new G4PVPlacement(0, G4ThreeVector(405.5*mm, -200.0, -159.0), Brcks10, "Brck10", detectorLV, false, 0);
    G4LogicalVolume * Brcks11=
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV11");
    new G4PVPlacement(0, G4ThreeVector(-1, -200.0, -159.0), Brcks11, "Brck11", detectorLV, false, 0);
    G4LogicalVolume * Brcks12 =
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV12");
    new G4PVPlacement(0, G4ThreeVector(-407.5, -200.0, -159.0), Brcks12,    "Brck12", detectorLV, false, 0);
    G4LogicalVolume * Brcks13=
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV13");
    new G4PVPlacement(0, G4ThreeVector(-814.0, -200.0, -159.0), Brcks13, "Brck13", detectorLV, false, 0);
    G4LogicalVolume * Brcks14 =
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV14");
    new G4PVPlacement(0, G4ThreeVector(812.0, -200.0, -159.0), Brcks14,    "Brck14", detectorLV, false, 0);
    
    ////Rectangle Shapes Right of Square brick
    G4Box * Rectngl =
    new G4Box("RctnglS", br_x, br_y, rct_l);
    G4LogicalVolume * Rctngl1=
    new G4LogicalVolume(Rectngl, fPb, "BrcksLV15");
    new G4PVPlacement(0, G4ThreeVector(407.0, -200.0, 451.0), Rctngl1, "Rectngl1", detectorLV, false, 0);
    G4LogicalVolume * Rctngl2 =
    new G4LogicalVolume(Rectngl, fPb, "BrcksLV16");
    new G4PVPlacement(0, G4ThreeVector(610.5, -200.0, 451.0), Rctngl2,    "Rectngl2", detectorLV, false, 0);
    G4LogicalVolume * Rctngl3=
    new G4LogicalVolume(Rectngl, fPb, "BrcksLV17");
    new G4PVPlacement(0, G4ThreeVector(814.0, -200.0, 451.0), Rctngl3, "Rectngl3", detectorLV, false, 0);
    ///Left Small Rectangle
    G4Box * Rectng2 =
    new G4Box("RctnglS2", br_x, br_y, bx_l);
    G4LogicalVolume * Rctngl4 =
    new G4LogicalVolume(Rectng2, fPb, "BrcksLV18");
    new G4PVPlacement(0, G4ThreeVector(-814.0, -200.0, 654.5), Rctngl4,    "Rectngl4", detectorLV, false, 0);
    ///Top Row of Bricks
    
    
    G4LogicalVolume * Brcks15=
    new G4LogicalVolume(Brcks, fPb, "BrcksLV19");
    new G4PVPlacement(0, G4ThreeVector(0.0,-200.0*mm,-464.5), Brcks15, "Brck15", detectorLV, false, 0);
    
    G4LogicalVolume * Brcks16 =
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV20");
    new G4PVPlacement(0, G4ThreeVector(305.0, -200.0, -362.5), Brcks16, "Brck16", detectorLV, false, 0);
    G4LogicalVolume * Brcks17=
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV21");
    new G4PVPlacement(0, G4ThreeVector(-305.0, -200.0, -362.5), Brcks17, "Brck17", detectorLV, false, 0);
    G4LogicalVolume * Brcks18 =
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV22");
    new G4PVPlacement(0, G4ThreeVector(305.0, -200.0, -566.0), Brcks18,    "Brck18", detectorLV, false, 0);
    G4LogicalVolume * Brcks19=
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV23");
    new G4PVPlacement(0, G4ThreeVector(-305.0, -200.0, -566.0), Brcks19, "Brck19", detectorLV, false, 0);
    G4LogicalVolume * Brcks20 =
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV24");
    new G4PVPlacement(0, G4ThreeVector(711.5, -200.0, -362.5), Brcks20,    "Brck20", detectorLV, false, 0);
    G4LogicalVolume * Brcks21=
    new G4LogicalVolume(BrcksSW, fPb, "BrcksLV25");
    new G4PVPlacement(0, G4ThreeVector(-711.5, -200.0, -362.5), Brcks21, "Brck21",  detectorLV, false, 0);
    
    
    return detectorLV;*/
}


G4LogicalVolume*GeometryConstruction::PentagonConstruct()
{
   /* const G4double br_x = 10.16*cm;
    const G4double br_y = 5.08*cm;
    const G4double br_z = 20.32*cm;
    
    G4RotationMatrix * rotY = new G4RotationMatrix();
    rotY->rotateY(60*deg);
    G4RotationMatrix * rotY2 = new G4RotationMatrix();
    rotY2->rotateY(120*deg);
    
    G4Box * detectorBox =
    new G4Box("detectorBox", 10000.0, 10000.0, 10000.0);
    G4LogicalVolume * detectorLV =
    new G4LogicalVolume(detectorBox, fAir, "detLV");
    // G4PVPlacement * detectorPhys =
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), detectorLV,
                      "detPhys", 0, false, 0);
    G4Box * BrcksP =
    new G4Box("Brcks2", br_z, br_y, br_x);
    
    ///Three Bottom Bricks
    G4LogicalVolume * Brcks22 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV26");
    new G4PVPlacement(0, G4ThreeVector(0.0, -200.0, 850.0), Brcks22, "Brck22", detectorLV, false, 0);
    G4LogicalVolume * Brcks23=
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV27");
    new G4PVPlacement(0, G4ThreeVector(406.5, -200.0, 850.0), Brcks23, "Brck23", detectorLV, false, 0);
    G4LogicalVolume * Brcks24 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV28");
    new G4PVPlacement(0, G4ThreeVector(-406.5, -200.0, 850.0), Brcks24,    "Brck24", detectorLV, false, 0);
    /////Three Right Diagonal Bricks
    G4LogicalVolume * Brcks25 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV29");
    new G4PVPlacement(rotY, G4ThreeVector(782.0, -200.0, 179), Brcks25, "Brck25", detectorLV, false, 0);
    G4LogicalVolume * Brcks26 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV30");
    new G4PVPlacement(rotY, G4ThreeVector(578.5, -200.0, -175.0), Brcks26,    "Brck26", detectorLV, false, 0);
    G4LogicalVolume * Brcks27 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV31");
    new G4PVPlacement(rotY, G4ThreeVector(374.0, -200.0, -529), Brcks27, "Brck27",  detectorLV, false, 0);
    //////Three Left Diagonal Bricks
    G4LogicalVolume * Brcks28 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV32");
    new G4PVPlacement(rotY2, G4ThreeVector(-782.0, -200.0, 179), Brcks28, "Brck28", detectorLV, false, 0);
    G4LogicalVolume * Brcks29 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV33");
    new G4PVPlacement(rotY2, G4ThreeVector(-578.5, -200.0, -175.0), Brcks29, "Brck29", detectorLV, false, 0);
    G4LogicalVolume * Brcks30 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV34");
    new G4PVPlacement(rotY2, G4ThreeVector(-374.0, -200.0, -529), Brcks30,    "Brck30", detectorLV, false, 0);
    //////Two Solo Diagonal and Solo Top Brick
    G4LogicalVolume * Brcks31 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV35");
    new G4PVPlacement(0, G4ThreeVector(0.0, -200.0, -519.0), Brcks31, "Brck31", detectorLV, false, 0);
    G4LogicalVolume * Brcks32 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV36");
    new G4PVPlacement(rotY2, G4ThreeVector(600.0, -200.0, 521), Brcks32,    "Brck32", detectorLV, false, 0);
    G4LogicalVolume * Brcks33 =
    new G4LogicalVolume(BrcksP, fPb, "BrcksLV37");
    new G4PVPlacement(rotY, G4ThreeVector(-600.0, -200.0, 521), Brcks33, "Brck33",  detectorLV, false, 0);
    
    return detectorLV;*/
}

void GeometryConstruction::SetReadFile( const G4String& File )
{
    fReadFile=File;
    fWritingChoice=0;
}
