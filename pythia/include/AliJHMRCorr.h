#ifndef ALIJHMRCORR_H
#define ALIJHMRCORR_H

#include "AliJBaseTrack.h"
#include "AliJHMRConst.h"
#include "AliJHMRHist.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <Pythia8/Pythia.h>

#include <TString.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH2D.h>

#include <vector>

using namespace std;
using namespace Pythia8;

class AliJHMRCorr {

public:

    AliJHMRCorr(AliJHMRHist *inhistos, detector det, bool isFullSim) :
        histos(inhistos) {
            fPhotonEfficiency = new TF1("fPhotonEfficiency", "TMath::Exp(-3.20093/x)"); // Parameters from fit to efficiency (PhotonEfficiency.C)
            fPhotonAcceptanceEfficiency = new TF1("fPhotonAcceptanceEfficiency", accFunc[det].Data()); // Parameters from fit (CheckMissingPionsRatio.C)
            std::cout << "Correlations are calculated in " << detEta[det][0]+etacut << " < eta < " << detEta[det][1]-etacut << std::endl;
            std::cout << "(etacut=" << etacut << ")" << std::endl;
            std::cout << "Using acceptance function " << accFunc[det].Data() << std::endl;
            fRand = new TRandom3();
            fIsFullSim = isFullSim;
        }

    virtual ~AliJHMRCorr(){ }

    bool IsTrackerAcceptance(double eta, double etaRange=etaTrackerRange);
    bool IsDetAcceptance(double eta, detector labelDet);
    int GetBin(double arr[], int nArr, double val);
    double GetDeltaPhi(double phiTrigg, double phiAssoc);
    double PhotonEnergySmearing(double px, double py, double pz);
    void SmearEnergies(TClonesArray * arrParticles);
    bool IsPhotonRemoved(double ePhoton);
    bool IsPhotonRemoved(TClonesArray *arrPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
    double GetAsymmetry(TClonesArray *arrPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
    AliJBaseTrack GetPhotonSumVector(TClonesArray *arrayPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
    int GetLeadingTriggerIndex(TClonesArray *arrPi0, bool bUseSim);
    int GetLargerTrigg(TClonesArray *arrPi0Peak, std::vector<int> listTriggPeak, TClonesArray *arrPi0Side, std::vector<int> listTriggSide);

    void DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bUseLeading, bool bUseWeight);
    void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bUseLeading, bool bUseWeightTrigg, bool bUseWeightAssoc);
    void ConstructTrueCorrComponents(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, bool bUseWeight);
    int ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, detector idet, bool bMass);
    void GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bMass, bool bUseLeading);

    bool IsMassWindow(double mass);
    bool IsMassWindow(double mass, int ibin, bool isTriggBin);
    bool IsSideband(double mass);

    void FillRealTriggers(TClonesArray *arrRealPi0, std::vector<int>& listTrigg);
    void FillPionMasses(TClonesArray *arrPhoton, int binsWithTriggPeak[NTRIGGBINS], int binsWithTriggSide[NTRIGGBINS], detector idet);
    void FillPionMassesLeading(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, std::vector<int> listTrigg, TH1D *hMass, detector idet);
    void FillPionMassesTrue(TClonesArray *arrPi0, int binsWithTriggReal[NTRIGGBINS], detector idet);
    void FillAsymmetry(TClonesArray *arrPhoton, detector idet);

protected:

    bool fIsFullSim; // choose if data is pythia monte carlo data or from geant simulation

    TF1 *fPhotonEfficiency;
    TF1 *fPhotonAcceptanceEfficiency;
    TRandom3 *fRand;
    AliJHMRHist *histos;
};

#endif
