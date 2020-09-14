#include "SteppingAction.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"
#include "RootIO.hh" 
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction()
: G4UserSteppingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step * theStep)
{
    
    // If you want to do anything for every step, do it here
    
    /*G4StepStatus stepStatus = theStep->GetPostStepPoint()->GetStepStatus();
    G4bool transmit = (stepStatus==fGeomBoundary || stepStatus==fWorldBoundary);
    if (transmit) {return;}
    //secondaries
    G4ParticleDefinition* particle = theStep->GetTrack()->GetDefinition();
    const std::vector<const G4Track*>* secondary  = theStep->GetSecondaryInCurrentStep();
    for (size_t lp=0; lp<(*secondary).size(); lp++) {
        particle = (*secondary)[lp]->GetDefinition();
        G4String name   = particle->GetParticleName();
        G4String type   = particle->GetParticleType();
        G4double energy = (*secondary)[lp]->GetKineticEnergy();
        G4ThreeVector pos = (*secondary)[lp]->GetPosition();
        G4double weight = (*secondary)[lp]->GetWeight();
        
        //G4int ih = 0;
        //if (particle == G4Gamma::Gamma())       ih = 2;
        //else if (particle == G4Neutron::Neutron())   ih = 3;
        //else if (particle == G4Proton::Proton())     ih = 4;
        //else if (particle == G4Deuteron::Deuteron()) ih = 5;
        //else if (particle == G4Alpha::Alpha())       ih = 6;
        //else if (type == "nucleus")                  ih = 7;
        //else if (type == "meson")                    ih = 8;
        //else if (type == "baryon")                   ih = 9;
        
        if (particle == G4Neutron::Neutron()){
            RootIO::GetInstance()->FillNeutronStuff(energy,pos.x(),pos.y(),pos.z(),weight);
        }
     
    }*/
    
}
