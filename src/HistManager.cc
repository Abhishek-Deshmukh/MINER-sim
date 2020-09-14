#include "HistManager.hh"
#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
#include "Cintex/Cintex.h"
#endif


HistManager::HistManager(TFile* myFile){
    aFile =  myFile;
}


HistManager::~HistManager(){
    
}


void HistManager::writeHists(){
    
    aFile->cd();
    std::map<std::string,TH1F*>::const_iterator mapit1;
    std::map<std::string,TH2F*>::const_iterator mapit2;
    std::map<std::string,TProfile*>::const_iterator mapit3;
    for (mapit1 = the1DMap.begin(); mapit1 != the1DMap.end(); ++mapit1){
        (*mapit1).second->Write();
    }
    for (mapit2 = the2DMap.begin(); mapit2 != the2DMap.end(); ++mapit2){
        (*mapit2).second->Write();
    }
    for (mapit3 = theProfMap.begin(); mapit3 != theProfMap.end(); ++mapit3){
        (*mapit3).second->Write();
    }
    aFile->cd();
    aFile->Close();
    
}


void HistManager::fill1DHist(float x, std::string name, std::string title,
                             int bins, float xmin, float xmax, float weight, std::string folder){
    
    std::map<std::string,TH1F*>::iterator it;
    it = the1DMap.find(name);
    if (it == the1DMap.end()){
        if (aFile->cd(folder.c_str())){ aFile->cd(folder.c_str());}
        else {aFile->mkdir(folder.c_str()); aFile->cd(folder.c_str());}
        the1DMap[name] = new TH1F(name.c_str(),title.c_str(),bins,xmin,xmax);
        the1DMap[name]->Sumw2();
        aFile->cd();
    }
    
    the1DMap[name]->Fill(x,weight);
    
}

void HistManager::fill1DHistUnevenBins(float x, std::string name, std::string title,
                                       int bins, float *binEdges, float weight, std::string folder){
    std::map<std::string,TH1F*>::iterator it;
    it = the1DMap.find(name);
    if (it == the1DMap.end()){
        if (aFile->cd(folder.c_str())){ aFile->cd(folder.c_str());}
        else {aFile->mkdir(folder.c_str()); aFile->cd(folder.c_str());}
        the1DMap[name] = new TH1F(name.c_str(),title.c_str(),bins,binEdges);
        the1DMap[name]->Sumw2();
        aFile->cd();
    }
    
    
    the1DMap[name]->Fill(x,weight);
    
}




void HistManager::fill2DHist(float x, float y, std::string name, std::string title,
                             int binsx, float xmin, float xmax,
                             int binsy, float ymin, float ymax, float weight, std::string folder){
    
    std::map<std::string,TH2F*>::iterator it;
    it = the2DMap.find(name);
    if (it == the2DMap.end()){
        if (aFile->cd(folder.c_str())){ aFile->cd(folder.c_str());}
        else {aFile->mkdir(folder.c_str()); aFile->cd(folder.c_str());}
        the2DMap[name] = new TH2F(name.c_str(),title.c_str(),binsx,xmin,xmax,binsy,ymin,ymax);
        the2DMap[name]->Sumw2();
        aFile->cd();
    }
    
    the2DMap[name]->Fill(x,y,weight);
    
}

void HistManager::fill2DHistUnevenBins(float x, float y, std::string name, std::string title,
                                       int binsx, float *binEdgesx,
                                       int binsy, float *binEdgesy, float weight,  std::string folder){
    
    std::map<std::string,TH2F*>::iterator it;
    it = the2DMap.find(name);
    if (it == the2DMap.end()){
        if (aFile->cd(folder.c_str())){ aFile->cd(folder.c_str());}
        else {aFile->mkdir(folder.c_str()); aFile->cd(folder.c_str());}
        the2DMap[name] = new TH2F(name.c_str(),title.c_str(),binsx,binEdgesx,binsy,binEdgesy);
        the2DMap[name]->Sumw2();
        aFile->cd();
    }
    
    
    the2DMap[name]->Fill(x,y,weight);
}


void HistManager::fill3DHist(float x, float y, float z, std::string name, std::string title,
                             int binsx, float xmin, float xmax,
                             int binsy, float ymin, float ymax,
                             int binsz, float zmin, float zmax,
                             float weight,  std::string folder){
    
    std::map<std::string,TH3F*>::iterator it;
    it = the3DMap.find(name);
    if (it == the3DMap.end()){
        if (aFile->cd(folder.c_str())){ aFile->cd(folder.c_str());}
        else {aFile->mkdir(folder.c_str()); aFile->cd(folder.c_str());}
        the3DMap[name] = new TH3F(name.c_str(),title.c_str(),binsx,xmin,xmax,binsy,ymin,ymax,binsz,zmin,zmax);
        the3DMap[name]->Sumw2();
        aFile->cd();
    }
    
    
    the3DMap[name]->Fill(x,y,z,weight);
}


void HistManager::fillProfile(float x, float y, std::string name, std::string title,
                              int binsx, float xmin, float xmax,
                              float ymin, float ymax, float weight,  std::string folder){
    
    std::map<std::string,TProfile*>::iterator it;
    it = theProfMap.find(name);
    if (it == theProfMap.end()){
        if (aFile->cd(folder.c_str())){ aFile->cd(folder.c_str());}
        else {aFile->mkdir(folder.c_str()); aFile->cd(folder.c_str());}
        theProfMap[name] = new TProfile(name.c_str(),title.c_str(),binsx,xmin,xmax,ymin,ymax);
        theProfMap[name]->Sumw2();
        aFile->cd();
    }
    
    theProfMap[name]->Fill(x,y,weight);
    
}
