#ifndef ALIJHMRHIST_H
#define ALIJHMRHIST_H

#include "AliJHMRConst.h"

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
        virtual ~AliJHMRHist(){ }    //destructor
        
        AliJHMRHist(const AliJHMRHist& obj);
        AliJHMRHist& operator=(const AliJHMRHist& obj);


        // create histograms 
        void CreateHistos();

    protected:

    //-------------------------
    //  Histogram parameters  |
    //-------------------------

    TH1D *hPionPt;
    TH1D *hPionPtFor;
    TH1D *hPionPtForDetected;
    TH1D *hPionPtMid;
    
    TH1D *hChargedHadronPt;
    TH1D *hChargedHadronPtFor;
    TH1D *hChargedHadronPtMid;

    TH1D *hPhotonPt;
    TH1D *hPhotonPtFor;
    TH1D *hPhotonPtMid;
    TH1D *hPhotonEnergyReal;
    TH1D *hPhotonEnergy;
    
    TH1D *hPionEta;
    TH1D *hChargedHadronEta;

    TH2D *hCorrMid[nTriggBins][nAssocBins];
    TH2D *hCorrFor[nTriggBins][nAssocBins];
    TH2D *hCorrChargedMid[nTriggBins][nAssocBins];
    TH2D *hCorrChargedFor[nTriggBins][nAssocBins];
    
    TH2D *hCorrMassMass[nTriggBins][nAssocBins];
    TH2D *hCorrMassSide[nTriggBins][nAssocBins];
    TH2D *hCorrSideMass[nTriggBins][nAssocBins];
    TH2D *hCorrSideSide[nTriggBins][nAssocBins];

    TH1D *hPi0MassTrigg[nTriggBins];
    TH1D *hPi0MassAssocPeak[nTriggBins][nAssocBins];
    TH1D *hPi0MassAssocSide[nTriggBins][nAssocBins];

    TDirectory *dirMasses;
    TDirectory *dirCorrMid;
    TDirectory *dirCorrFor;
    TDirectory *dirCorrChargedMid;
    TDirectory *dirCorrChargedFor;
    TDirectory *dirCorrMassMass;
    TDirectory *dirCorrMassSide;
    TDirectory *dirCorrSideMass;
    TDirectory *dirCorrSideSide;
};

#endif























