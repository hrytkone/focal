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
        void FillMathingInformation(TClonesArray * arrClusters, TClonesArray * arrPtMatchedClusters);

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
        TH2D *hCorrRejectMassMass[NTRIGGBINS][NASSOCBINS]; // Control histogram for rejected pairs in mass-mass situation
        TH2D *hCorrRejectSideSide[NTRIGGBINS][NASSOCBINS]; // Control histogram for rejected pairs in side-side situation
        TH2D *hCorrChargedMid[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrChargedFor[NTRIGGBINS][NASSOCBINS];

        TH2D *hCorrMassMass[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrMassSide[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideMass[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideSide[NTRIGGBINS][NASSOCBINS];

        // To check the different sources in the background
        TH2D *hCorrSideSideDecay[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideSideMix[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideSideNotDecay[NTRIGGBINS][NASSOCBINS];

        TH2D *hCorrMassMassMixed[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrMassSideMixed[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideMassMixed[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSideSideMixed[NTRIGGBINS][NASSOCBINS];

        TH2D *hCorrSignalSignal[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrSignalBg[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrBgSignal[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrBgBg[NTRIGGBINS][NASSOCBINS];

        // To check the different sources in the background
        TH2D *hCorrBgBgDecay[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrBgBgMix[NTRIGGBINS][NASSOCBINS];
        TH2D *hCorrBgBgNotDecay[NTRIGGBINS][NASSOCBINS];

        TH1D *hMassTrue;
        TH1D *hMassDecay;
        TH1D *hMassMix;
        TH1D *hMassNotDecay;

        TH1D *hPi0MassTrigg[NTRIGGBINS];
        TH1D *hPi0MassAssocPeak[NTRIGGBINS][NASSOCBINS];
        TH1D *hPi0MassAssocSide[NTRIGGBINS][NASSOCBINS];

        TH1D *hEnergyAsymTrue;
        TH1D *hEnergyAsymRec;

        TH2D *hMassAsymTrue[NASSOCBINS];
        TH2D *hMassAsymDecay[NASSOCBINS];
        TH2D *hMassAsymMix[NASSOCBINS];
        TH2D *hMassAsymNotDecay[NASSOCBINS];
        TH2D *hMassOpeningAngleTrue[NASSOCBINS];
        TH2D *hMassOpeningAngleDecay[NASSOCBINS];
        TH2D *hMassOpeningAngleMix[NASSOCBINS];
        TH2D *hMassOpeningAngleNotDecay[NASSOCBINS];

        // To check how well matching particles to clusters work
        TH2D *hPtMatched;
        TH2D *hEtaMatched;
        TH2D *hPhiMatched;

        TH1D *hEnergyMassBgTrigg[NTRIGGBINS];
        TH1D *hEnergySidebandTrigg[NTRIGGBINS];
        TH1D *hEnergyMassBgAssoc[NTRIGGBINS][NASSOCBINS];
        TH1D *hEnergySidebandAssoc[NTRIGGBINS][NASSOCBINS];

        TDirectory *dirMasses;
        TDirectory *dirMassComponents;
        TDirectory *dirCorrMid;
        TDirectory *dirCorrFor;
        TDirectory *dirCorrMeas;
        TDirectory *dirCorrReject;
        TDirectory *dirCorrChargedMid;
        TDirectory *dirCorrChargedFor;
        TDirectory *dirCorrMassMass;
        TDirectory *dirCorrMassSide;
        TDirectory *dirCorrSideMass;
        TDirectory *dirCorrSideSide;
        TDirectory *dirTrueComponents; // true components f_SS, f_SB, f_BS, f_BB
        TDirectory *dirBackgroundSources;
        TDirectory *dirBackgroundEnergies;
        TDirectory *dirAsym;
};

#endif
