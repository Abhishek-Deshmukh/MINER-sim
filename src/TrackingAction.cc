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
// $Id: TrackingAction.cc 81501 2014-05-31 13:51:51Z ldesorgh $
//
/// \file eventgenerator/exgps/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//

#include "TrackingAction.hh"
#include "TrackExtra.hh"

#include "G4Track.hh"
//#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "RootIO.hh"
#include "G4TrackingManager.hh"

TrackingAction::TrackingAction()
:G4UserTrackingAction()
{ }

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{

    if (track->GetParentID() == 0){
        RootIO::GetInstance()->SetIncomingE(track->GetKineticEnergy());
        RootIO::GetInstance()->AddTrack(track);
        //G4Track *mutableTrack = fpTrackingManager->GetTrack();
        //TrackExtra *info = new TrackExtra(track);
        //mutableTrack->SetUserInformation(info);
    }


    // These are tracks created by the importance sampling
    /*
     if (track->GetParentID() != 0){
     if (  track->GetCreatorProcess()->GetProcessName() == "ImportanceProcess"){
     RootIO::GetInstance()->AddTrack(track);
     }
     }
     */

}

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
    /*
     G4Track *mutableTrack = fpTrackingManager->GetTrack();
     if (track->GetParentID() == 0){
     G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
     if(secondaries)
     {
     TrackExtra* info = (TrackExtra*)(track->GetUserInformation());
     size_t nSeco = secondaries->size();
     if(nSeco>0)
     {
     for(size_t i=0;i<nSeco;i++)
     {
     TrackExtra* infoNew = new TrackExtra(info);
     (*secondaries)[i]->SetUserInformation(infoNew);
     }
     }
     }
     }
     if (track->GetParentID() != 0){
     if (  track->GetCreatorProcess()->GetProcessName() == "ImportanceProcess"){
     TrackExtra *info = new TrackExtra(track);
     mutableTrack->SetUserInformation(info);
     G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
     if(secondaries)
     {
     size_t nSeco = secondaries->size();
     if(nSeco>0)
     {
     for(size_t i=0;i<nSeco;i++)
     {
     TrackExtra* infoNew = new TrackExtra(info);
     (*secondaries)[i]->SetUserInformation(infoNew);
     }
     }
     }
     }
     else {
     G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
     if(secondaries)
     {
     TrackExtra* info = (TrackExtra*)(track->GetUserInformation());
     size_t nSeco = secondaries->size();
     if(nSeco>0)
     {
     for(size_t i=0;i<nSeco;i++)
     {
     TrackExtra* infoNew = new TrackExtra(info);
     (*secondaries)[i]->SetUserInformation(infoNew);
     }
     }
     }
     }
     }
     */
}
