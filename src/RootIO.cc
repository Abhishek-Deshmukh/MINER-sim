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
/// \file persistency/P01/src/RootIO.cc
/// \brief Implementation of the RootIO class
//
// $Id: RootIO.cc 82130 2014-06-11 09:26:44Z gcosmo $
//

#include <sstream>

#include "RootIO.hh"
#include "RootIOMessenger.hh"


//
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Threading.hh"
#include "TrackExtra.hh"


#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
#include "Cintex/Cintex.h"
#endif


static RootIO* instance = 0;

RootIO::RootIO():fileName("hits.root")
{

    fMessenger = new RootIOMessenger(this);

}

RootIO::~RootIO()
{}

void RootIO::Setup()
{
    //G4int thread = G4Threading::G4GetThreadId();
    //std::ostringstream temp;
    //if (thread >= 0){
    //  temp<<fileName<<"_"<<thread<<".root";
    //} else { temp<<fileName<<".root"; }
    theFile = new TFile(fileName,"RECREATE");
    theFile->cd();
    theTree  = new TTree("theTree", "go climb one");
    //sHits  = new TClonesArray("BaseHit");
    //theTree->Branch("sHits",&sHits, 6400, 0);
    hitC = 0;
    //sTracks  = new TClonesArray("BaseTrack");
    //theTree->Branch("sTracks",&sTracks, 6400, 0);
    trackC = 0;

    theTree->Branch("hitC",&hitC,"hitC/I");
    theTree->Branch("trackC",&trackC,"trackC/I");
    theTree->Branch("event",&event,"event/I");
    theTree->Branch("eInc",&eInc,"eInc/F");

    hitC = 0;
    trackC = 0;
    event = 1;
    eInc = 0.;

    hists = new HistManager(theFile);

}

void RootIO::SetFileName(G4String file)
{
    fileName = file;
}

G4String RootIO::GetFileName()
{
    return fileName;
}

RootIO* RootIO::GetInstance()
{
    if (instance == 0 )
    {
        instance = new RootIO();
    }
    return instance;
}

void RootIO::SetIncomingE(G4double en)
{

    eInc = en;

}

void RootIO::AddHits(MinerHitsCollection * zipHits, G4int detID)
{
    for (G4int i = 0; i < zipHits->entries(); i++)
    {
        /*MinerHit* hit = new ((*sHits)[hitC]) MinerHit;
        G4ThreeVector vec = (*zipHits)[i]->GetPos();
        G4ThreeVector mom = (*zipHits)[i]->GetMom();
        hit->SetPos(vec.x(),vec.y(),vec.z());
        hit->Setp3(mom.x(),mom.y(),mom.z());
        hit->SetEkin((*zipHits)[i]->GetParticleEnergy());
        hit->SetEdep((*zipHits)[i]->GetEdep());
        hit->SetPid((*zipHits)[i]->GetPDGID());
        hit->SetDetID(detID);
        hit->SetTime((*zipHits)[i]->GetTime());
        hit->SetWeight((*zipHits)[i]->GetWeight());
        hit->SetPreProcess((*zipHits)[i]->GetPreProcess());
        hit->SetPostProcess((*zipHits)[i]->GetPostProcess());
        hit->SetInc(0);*/


        /*G4double edep = zipHits->GetTotalEnergyDeposit();



        MinerHit* newHit = new MinerHit();


        G4int postproc = 0;



        newHit->SetMom(zipHits[i]->GetTrack()->GetMomentum());
        newHit->SetTrackID  (zipHits[i]->GetTrack()->GetTrackID());
        newHit->SetPDGID (zipHits[i]->GetTrack()->GetDefinition()->GetPDGEncoding());
        newHit->SetTime(zipHits[i]->GetPostStepPoint()->GetGlobalTime());
        newHit->SetEdep(zipHits[i]->GetTotalEnergyDeposit());
        newHit->SetPos (zipHits[i]->GetPostStepPoint()->GetPosition());
        newHit->SetParticleEnergy(zipHits[i]->GetPreStepPoint()->GetKineticEnergy());
        newHit->SetWeight(zipHits[i]->GetTrack()->GetWeight());
        newHit->SetPreProcess(postproc);
        newHit->SetPostProcess(postproc);

        */

   }
}

void RootIO::AddTrack(const G4Track*  trk)
{

    //G4double ekin           = trk->GetKineticEnergy();
    //G4double weight         = trk->GetWeight();



    G4ThreeVector vertex = trk->GetVertexPosition ();
    //G4double x = vertex.x(), y = vertex.y(), z = vertex.z();

    G4ThreeVector mom = trk->GetMomentum();
    G4double px = mom.x(), py = mom.y(), pz = mom.z(), pt = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

    hists->fill1DHist(px/pt,"Primary_px","",100,-1.1,1.1,1.,"");
    hists->fill1DHist(py/pt,"Primary_py","",100,-1.1,1.1,1.,"");
    hists->fill1DHist(pz/pt,"Primary_pz","",100,-1.1,1.1,1.,"");

    std::cout<<"px is "<<px<<" py is "<<py<<" pz is "<<pz<< " ptot is "<< pt <<std::endl;

   // BaseTrack* tr = new ((*sTracks)[trackC]) BaseTrack;
    //tr->Setp4(px,py,pz,trk->GetTotalEnergy());
    //tr->SetPos(x,y,z);
    //tr->SetPid(trk->GetDynamicParticle()->GetPDGcode());
   // tr->SetTime(trk->GetGlobalTime());
    //tr->SetInc(trk->GetParentID());
    //trackC++;

}

void RootIO::FillMonitoring(MinerHitsCollection * zipHits, G4int detID)
{

    float edepo = 0;

    std::ostringstream detNumber_ss;

    detNumber_ss << detID;
    G4String detNumber = detNumber_ss.str();

    for (G4int i = 0; i < zipHits->entries(); i++)
    {
        G4ThreeVector vec = (*zipHits)[i]->GetPos();
        G4ThreeVector mom = (*zipHits)[i]->GetMom();

        // Recomendation: use std::to_string() if you can require C++11 (which Geant 10.2 does)

        std::ostringstream pidNumber_ss;


        pidNumber_ss << (*zipHits)[i]->GetPDGID();

        G4String pidNumber = pidNumber_ss.str();

        edepo = edepo + (*zipHits)[i]->GetEdep();

        hists->fill1DHist(edepo,"Det"+detNumber+"_Ekin_PID"+pidNumber,"",10000,0,2,1,"Det"+detNumber+"_Monitoring");

        hists->fill2DHist(vec.y(),vec.z(),"Det"+detNumber+"_pos_PID"+pidNumber,"",120,-600,600,120,-1900,-700,(*zipHits)[i]->GetWeight(),"Det"+detNumber+"_Monitoring");


        //hists->fill3DHist(vec.x(),vec.y(), vec.z(), "Hit_Location_"+detNumber, "Hit_Location",
        //                 100, -1000, 1000,
        //                 100, -1000, 1000,
        //                 100, -1000, 1000,
        //                 1, "Det"+detNumber);


      //if(pidNumber=="12"){
      //      std::cout<<pidNumber+"Energy is"+edepo<<std::endl;
      //      std::cout<<detNumber<<std::endl;
      //  }


    }

     hists->fill1DHist(edepo*1000,"Energy_Deposited"+detNumber,"Energy Deposited",4096,0,5250,1,"Det"+detNumber);


    /*for(int j = 0; j < zipHits->entries(); j++)
    {


        edepo = edepo + (*zipHits)[j]->GetEdep();

    }

    std::ofstream outFile("edep.txt");
    G4double a = edepo;
    outFile << a << G4endl;


    hists->fill1DHist(edepo,"Energy Deposited","",500,0,2,1,"Det");

 */




}


void RootIO::FillNeutronStuff(G4double energy, G4double x, G4double y, G4double z, G4double weight)
{

    hists->fill1DHist(x,"Neutron_x","Neutron_x",4000,-1000,3000,weight,"");
    hists->fill1DHist(energy,"Neutron_E","Neutron_E",500,0,10,1.,"");
    hists->fill2DHist(y,z,"Neutron_yz","",120,-600,600,120,-1900,-700,weight,"");
    hists->fill2DHist(x,z,"Neutron_xz","",120,-600,600,400,-1000,3000,weight,"");

}

void RootIO::Write()
{

    event++;
    if (hitC > 0){
        theTree->Fill();
    }
    //sHits->Clear();
    //sTracks->Clear();
    hitC = 0;
    trackC = 0;
}

void RootIO::Close()
{
    theFile->cd();
    theFile->Write();
    theFile->Close();
}
