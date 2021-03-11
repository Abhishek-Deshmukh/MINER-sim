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
// $Id: PrimaryGeneratorAction.cc 67272 2013-02-13 12:05:47Z ihrivnac $
//

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include <iomanip>
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
#include "G4SystemOfUnits.hh"
//#include "radsource.h"
//#include "cpp_api.h"


PrimaryGeneratorAction::PrimaryGeneratorAction(const char *inputfile)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fParticleSource(0),
  SourceZPosition(0)
{
  fParticleGun    = new G4ParticleGun();
  fParticleSource = new G4GeneralParticleSource();
  std::ifstream inputFile;
  inputFile.open(inputfile,std::ios::in);
  char buffer[1000];

  if (inputFile.fail()) {
    if( *inputfile !=0)  //....only complain if a filename was given
      G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << G4endl;
    InputState=-1;
  }else{
    std::string setupString("");
    while ( !inputFile.getline(buffer,1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    //CRYSetup *setup=new CRYSetup(setupString,"cosmic_data");
    //CRYSetup *setup=new CRYSetup(setupString,"/home/samir/cry_v1.7/data/cosmics_0.data");
    CRYSetup *setup=new CRYSetup(setupString,"/home/mercury/Downloads/CRY/cry_v1.7/data");

    gen = new CRYGenerator(setup);

    // set random number generator
    // RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
    // setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState=0;
  }
  // create a vector to store the CRY particle properties
  vect=new std::vector<CRYParticle*>;

  // Create the table containing all particle names
  particleTable = G4ParticleTable::GetParticleTable();
  fMessenger = new PrimaryGeneratorMessenger(this);
}

void PrimaryGeneratorAction::UpdateCRY(std::string* MessInput)
{
  CRYSetup *setup=new CRYSetup(*MessInput,"/home/mercury/Downloads/CRY/cry_v1.7/data");
  gen = new CRYGenerator(setup);
  // set random number generator
  // RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
  // setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
  InputState=0;

}


void PrimaryGeneratorAction::UpdateCRYFile(G4String newValue)
{
  G4cout<<"Input file is=> "<<newValue<<G4endl;
  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(newValue,std::ios::in);
  char buffer[1000];

  if (inputFile.fail()) {
    G4cout << "Failed to open input file " << newValue << G4endl;
    G4cout << "Make sure to define the cry library on the command line" << G4endl;
    InputState=-1;
  }else{
    std::string setupString("");
    while ( !inputFile.getline(buffer,1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup=new CRYSetup(setupString,"/home/mercury/Downloads/CRY/cry_v1.7/data");
    gen = new CRYGenerator(setup);

    // set random number generator
    // RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
    // setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState=0;
  }
}
void PrimaryGeneratorAction::SetSourceZPosition(G4double dist)
{
  SourceZPosition = dist;
}
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fParticleSource;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

if (InputState != 0) {
    G4String* str = new G4String("CRY library was not successfully initialized");
    //G4Exception(*str);
    G4Exception("CryPrimaryGeneratorAction", "1",
                RunMustBeAborted, *str);
  }
  G4String particleName;
  vect->clear();
  gen->genEvent(vect);

  //....debug output
  // G4cout << "\nEvent=" << anEvent->GetEventID() << " "
  //        << "CRY generated nparticles=" << vect->size()
  //        << G4endl;
  for ( unsigned j=0; j<vect->size(); j++) {
    particleName=CRYUtils::partName((*vect)[j]->id());

    //if( particleName == "electron" ) continue;



    //....debug output
    // G4cout << "  "          << particleName << " "
    // 	   << "charge="      << (*vect)[j]->charge() << " "
    // 	   << "energy (MeV)=" << (*vect)[j]->ke()*MeV << " "
    // 	   << "pos (m)"
    // 	   << G4ThreeVector((*vect)[j]->x(), (*vect)[j]->y(), (*vect)[j]->z())
    // 	   << " " << "direction cosines "
    // 	   << G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w())
    // 	   << " " << G4endl;

    fParticleGun->SetParticleDefinition(particleTable->FindParticle((*vect)[j]->PDGid()));
    fParticleGun->SetParticleEnergy((*vect)[j]->ke()*MeV);
    fParticleGun->SetParticlePosition(G4ThreeVector((*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->z()*m + SourceZPosition));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w()));
    fParticleGun->SetParticleTime((*vect)[j]->t());
    fParticleGun->GeneratePrimaryVertex(anEvent);
    delete (*vect)[j];
  }
  //fParticleSource->GeneratePrimaryVertex(anEvent);
}
