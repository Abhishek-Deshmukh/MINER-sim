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
// $Id: PrimaryGeneratorAction.hh 71200 2013-06-12 08:19:59Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4ParticleTable.hh"
// #include "CRYSetup.h"
// #include "CRYGenerator.h"
// #include "CRYParticle.h"
// #include "CRYUtils.h"
#include "vector"
#include "RNGWrapper.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
//#include "PrimaryGeneratorMessenger.hh"

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class CRYGenerator;
class CRYParticle;
class PrimaryGeneratorMessenger;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(const char *inputfile);
  ~PrimaryGeneratorAction();
    
  virtual void GeneratePrimaries(G4Event*);
  void UpdateCRY(std::string* MessInput);
  void UpdateCRYFile(G4String newValue);
  void SetSourceZPosition(G4double);
private:
  G4ParticleGun* fParticleGun;
  G4GeneralParticleSource* fParticleSource;
  std::vector<CRYParticle*> *vect; // vector of generated particles
  G4ParticleTable* particleTable;
  CRYGenerator* gen;
  G4int InputState;
  G4double SourceZPosition;
  PrimaryGeneratorMessenger *fMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
