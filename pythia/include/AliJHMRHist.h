#ifndef ALIJHMRHIST_H
#define ALIJHMRHIST_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>

class AliJHMRHist {

    public:
        AliJHMRHist(); //constructor
        virtual ~AliJHMRHist(){ delete fhistoList; }    //destructor
        
        AliJHMRHist(const AliJHistos& obj);
        AliJHMRHist& operator=(const AliJHistos& obj);
        
        void ScaleInvar(){
            scaleNotEquidistantHisto( hChargedPtpp);
            scaleNotEquidistantHisto( hChargedPtAA);
        }

        void scaleNotEquidistantHisto(TH1D *hid, double sc=1){
            for(int i=1;i<= hid->GetNbinsX();i++){
                hid->SetBinContent(i,hid->GetBinContent(i)*sc/hid->GetBinWidth(i));
                hid->SetBinError(i,hid->GetBinError(i)*sc/hid->GetBinWidth(i));
            }   
        }

        // create histograms 
        void CreatePi0MassHistos();
        void CreateAzimuthCorrHistos();
        void CreateToyHistos();
        void CreateIAAMoons();
        void CreateXEHistos();

        void CreateEventTrackHistos();
        void CreateJetHistos();

        void CreatePairPtCosThetaStar();
        void CreatePtCorrHistos();
        void CreateRunByRunHistos(int runID, int runcounter);

        void ReadInclusiveHistos(const char *inclusFileName);
        
        TList *GetHistoList() { return fhistoList; }

        void UseDirectory(bool b) { fUseDirectory=b; }
        bool UseDrectory(){ return fUseDirectory; }

        TDirectory * MakeDirectory(TString name){
            JumpToDefalutDirectory();
            TDirectory * dir = gDirectory->GetDirectory(name);
            if( !dir ) dir = gDirectory->mkdir(name);
            dir->cd();
            return dir;
        }
        TDirectory * JumpToDefalutDirectory(){
            fTopDirectory->cd();
            return gDirectory;
        }

protected:

    //-------------------------
    //  Histogram parameters  |
    //-------------------------
    const double etaBinWidth = 0.025;
    const double phiBinWidth = 0.025;

    const double etaTrackerRange = 0.9;
    const double etaFocalMin = 3.2;
    const double etaFocalMax = 5.8;
    const double etaFocalRange = etaFocalMax - etaFocalMin;
    const int nEtaBinTracker = 2.*int(etaTrackerRange/etaBinWidth) + 1;
    const int nEtaBinFocal = int(etaFocalRange/etaBinWidth) + 1;

    const double deltaPhiMin = -TMath::Pi()/2.0;
    const double deltaPhiMax = 3.0*TMath::Pi()/2.0;
    const int nPhiBin = int((deltaPhiMax-deltaPhiMin)/phiBinWidth) + 1;

    const int nIncPtBin = 150;
    double logBinsX[nIncPtBin+1], limMin = 0.1, limMax = 100;
    const double logBW = (log(limMax) - log(limMin))/nIncPtBin;

    const int nIncEtaBin = 150;
    const double incEtaRange = 20.0;

    const int nPhotonEnergyBin = 150;
    double limPhotonEnergyMin = 0., limPhotonEnergyMax = 1500.;

};

#endif























