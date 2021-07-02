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
        AliJHMRHist() { } //constructor
        virtual ~AliJHMRHist(){ }    //destructor
        
        AliJHMRHist(const AliJHMRHist& obj);
        AliJHMRHist& operator=(const AliJHMRHist& obj);

        // create histograms 
        void CreateHistos(TFile *output);
        
        TH1D *hCounter;
        TH1D *hRealTriggCounter;

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

        TH2D *hCorrMid[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrFor[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrChargedMid[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrChargedFor[NTRIGGBINS][NASSOCBINS];
        
        TH2D *hCorrMassMass[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrMassSide[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideMass[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideSide[NTRIGGBINS][NASSOCBINS];

        TH1D *hPi0MassTrigg[NTRIGGBINS];
        TH1D *hPi0MassAssocPeak[NTRIGGBINS][NASSOCBINS];
        TH1D *hPi0MassAssocSide[NTRIGGBINS][NASSOCBINS];

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























