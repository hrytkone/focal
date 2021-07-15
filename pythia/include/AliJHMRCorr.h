#ifndef ALIJHMRCORR_H
#define ALIJHMRCORR_H

#include "AliJHMRConst.h"
#include "AliJHMRHist.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <Pythia8/Pythia.h>

#include <TString.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH2D.h>

#include <vector>

using namespace std;
using namespace Pythia8;

class AliJHMRCorr {
  
public:
    
    AliJHMRCorr(AliJHMRHist *inhistos) : 
        histos(inhistos) { 
            fPhotonEfficiency = new TF1("fPhotonEfficiency", "TMath::Exp(-3.20093/x)"); // Parameters from fit to efficiency (PhotonEfficiency.C)
            fPhotonAcceptanceEfficiency = new TF1("fPhotonAcceptanceEfficiency", "TMath::Exp(-0.117082/(x + 0.0832931))"); // Parameters from fit (CheckMissingPionsRatio.C)
            fRand = new TRandom3();
        }

    virtual ~AliJHMRCorr(){ }

    bool IsTrackerAcceptance(double eta, double etaRange=etaTrackerRange);
    bool IsFocalAcceptance(double eta, double etaMin=etaFocalMin, double etaMax=etaFocalMax);
    int GetBin(double arr[], int nArr, double val);
    double GetDeltaPhi(double phiTrigg, double phiAssoc);
    double PhotonEnergySmearing(double px, double py, double pz);
    void SmearEnergies(TClonesArray * arrParticles);
    bool IsPhotonRemoved(double ePhoton);
    TLorentzVector GetPhotonSumVector(TClonesArray *arrayPhoton, int iPhoton1, int iPhoton2);
    int GetLeadingTriggerIndex(TClonesArray *arrPi0);
    int GetLargerTrigg(TClonesArray *arrPi0Peak, std::vector<int> listTriggPeak, TClonesArray *arrPi0Side, std::vector<int> listTriggSide);

    void DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bUseWeight);
    void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bUseWeightTrigg, bool bUseWeightAssoc);
    void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, bool bMass);
    void GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bUseLeading);

    bool IsMassWindow(double mass);
    bool IsSideband(double mass);

    void FillRealTriggers(TClonesArray *arrRealPi0, std::vector<int>& listTrigg);
    void FillPionMasses(TClonesArray *arrPhoton, int binsWithTriggPeak[NTRIGGBINS], int binsWithTriggSide[NTRIGGBINS]);

protected:

    TF1 *fPhotonEfficiency;
    TF1 *fPhotonAcceptanceEfficiency;
    TRandom3 *fRand;
    AliJHMRHist *histos;
};

#endif

