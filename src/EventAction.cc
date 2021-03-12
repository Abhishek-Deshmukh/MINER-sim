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

#include "EventAction.hh"

#include "TClonesArray.h"


#include "MinerHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"

#include "RootIO.hh"


EventAction::EventAction()
: G4UserEventAction()
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{}

void EventAction::EndOfEventAction(const G4Event* event)
{

    G4SDManager* fSDM = G4SDManager::GetSDMpointer();
    G4HCofThisEvent* HCofEvent = event->GetHCofThisEvent();

    // Get the Hits for the hybrid
    MinerHitsCollection*  HybridHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("Hybrid_hits")));

    // Get the Hits for the HV
    MinerHitsCollection*  HVHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("HV_hits")));

    // Get the Hits for the CsI
    MinerHitsCollection*  CsIHits = (MinerHitsCollection*)(HCofEvent->GetHC(fSDM->GetCollectionID("CsI_hits")));



    if (HybridHits->entries() > 0){
    RootIO::GetInstance()->FillMonitoring(HybridHits,1);
    }

    if (HVHits->entries() > 0){
    RootIO::GetInstance()->FillMonitoring(HVHits,2);
    }

    if (CsIHits->entries() > 0){
    RootIO::GetInstance()->FillMonitoring(CsIHits,3);
    }


    RootIO::GetInstance()->Write();


}
