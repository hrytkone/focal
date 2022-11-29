#ifndef ALIJHMRHIST_H
#define ALIJHMRHIST_H

#include "AliJHMRConst.h"
#include "AliJBaseTrack.h"

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
#include <TClonesArray.h>

class AliJHMRHist {

    public:
        AliJHMRHist() { } //constructor
        virtual ~AliJHMRHist(){ }    //destructor

        AliJHMRHist(const AliJHMRHist& obj);
        AliJHMRHist& operator=(const AliJHMRHist& obj);

        // create histograms
        void CreateHistos(TFile *output, detector labelDet, bool bUseLeading);
        void FillPtEta(particleType itype, TClonesArray * arrParticles);

        TH1D *hCounter;
        TH1D *hRealTriggCounter;

        TH1D *hPionPt;
        TH1D *hRecPionPt;
        TH1D *hChargedHadronPt;
        TH1D *hPhotonPt;

        TH1D *hPionEta;
        TH1D *hRecPionEta;
        TH1D *hPhotonEta;
        TH1D *hChargedHadronEta;

        TH2D *hCorrMid[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrFor[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrMeas[NTRIGGBINS][NASSOCBINS]; // This is not weighted
        TH2D *hCorrChargedMid[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrChargedFor[NTRIGGBINS][NASSOCBINS];

        TH2D *hCorrMassMass[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrMassSide[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideMass[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideSide[NTRIGGBINS][NASSOCBINS];

        TH2D *hCorrMassMassMixed[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrMassSideMixed[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideMassMixed[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideSideMixed[NTRIGGBINS][NASSOCBINS];

        TH2D *hCorrSignalSignal[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSignalBg[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrBgSignal[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrBgBg[NTRIGGBINS][NASSOCBINS];

        TH1D *hPi0MassTrigg[NTRIGGBINS];
        TH1D *hPi0MassAssocPeak[NTRIGGBINS][NASSOCBINS];
        TH1D *hPi0MassAssocSide[NTRIGGBINS][NASSOCBINS];

        TH1D *hEnergyAsymTrue;
        TH1D *hEnergyAsymRec;

        TDirectory *dirMasses;
        TDirectory *dirCorrMid;
        TDirectory *dirCorrFor;
        TDirectory *dirCorrMeas;
        TDirectory *dirCorrChargedMid;
        TDirectory *dirCorrChargedFor;
        TDirectory *dirCorrMassMass;
        TDirectory *dirCorrMassSide;
        TDirectory *dirCorrSideMass;
        TDirectory *dirCorrSideSide;
        TDirectory *dirTrueComponents; // true components f_SS, f_SB, f_BS, f_BB

};

#endif
