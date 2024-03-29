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

    AliJHMRCorr(AliJHMRHist *inhistos, detector det, bool bUseLeading, bool isFullSim) :
        histos(inhistos) {
            std::cout << "Correlations are calculated in " << detEta[det][0]+etacut << " < eta < " << detEta[det][1]-etacut << std::endl;
            std::cout << "(etacut=" << etacut << ")" << std::endl;

            fRand = new TRandom3();
            fUseLeading = bUseLeading;
            fIsFullSim = isFullSim;

            fLeta = -1.;
            fLphi = -1.;
            fLpt  = -1.;

            fX1 = -1.;
            fX2 = -1.;
        }

    virtual ~AliJHMRCorr(){ }

    bool IsTrackerAcceptance(double eta, double etaRange=etaTrackerRange);
    bool IsDetAcceptance(double eta, detector labelDet);
    int GetBin(double arr[], int nArr, double val);
    double GetDeltaPhi(double phiTrigg, double phiAssoc);
    double PhotonEnergySmearing(double px, double py, double pz);
    void SmearEnergies(TClonesArray * arrParticles);
    double GetAsymmetry(TClonesArray *arrPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
    double GetOpeningAngle(TClonesArray *arrPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
    AliJBaseTrack GetPhotonSumVector(TClonesArray *arrayPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
    int GetLeadingTriggerIndex(TClonesArray *arrPi0, bool bUseSim);
    int GetLargerTrigg(TClonesArray *arrPi0Peak, std::vector<int> listTriggPeak, TClonesArray *arrPi0Side, std::vector<int> listTriggSide);

    void DoCorrelations(TClonesArray *arrPi0, TClonesArray *arrPhoton, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bTrueCorr, bool bMassWindowTrigg, bool bUseWeight);
    void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bUseWeight);
    void ConstructTrueCorrComponents(TClonesArray *arrPi0, TClonesArray* arrPhoton, std::vector<int> listTrigg, std::vector<int> listAssoc, bool bUseWeight);
    int ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, detector idet, bool bMass);
    void GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bMass);

    bool IsMassWindow(double mass);
    bool IsMassWindow(double mass, int ibin, bool isTriggBin);
    bool IsSideband(double mass);
    bool IsLeadingTrigger(double mass, double pt);

    void FillRealTriggers(TClonesArray *arrRealPi0, std::vector<int>& listTrigg);
    void FillPionMasses(TClonesArray *arrPhoton, int binsWithTriggPeak[NTRIGGBINS], int binsWithTriggSide[NTRIGGBINS], detector idet);
    void FillPionMassesLeading(TClonesArray *arrPhoton, TClonesArray *arrPi0, std::vector<int> listTrigg, detector idet, bool isPeakTriggLarger);
    void FillPionMassesTrue(TClonesArray *arrPi0, int binsWithTriggReal[NTRIGGBINS], detector idet);
    void FillAsymmetry(TClonesArray *arrPhoton, detector idet);

    bool CheckAssocPhotonPair(int iTrigg, int iAssoc, bool bMassWindow);
    void ClearPhotonPairVector() { photonId.clear(); }
    void ClearSidebandPairVector() { sidebandId.clear(); }

    // For saving x1 and x2 distributions
    void SetX1X2(double x1, double x2) { fX1 = x1; fX2 = x2; }
    double GetX1() { return fX1; }
    double GetX2() { return fX2; }

protected:

    bool fUseLeading; // use leading particle correlation
    bool fIsFullSim;  // choose if data is pythia monte carlo data or from geant simulation

    double fLeta;     // These are to find out the leading trigger when filling mass histos for assoc
    double fLphi;     // These are to find out the leading trigger when filling mass histos for assoc
    double fLpt;      // These are to find out the leading trigger when filling mass histos for assoc

    double fX1;
    double fX2;

    TRandom3 *fRand;
    AliJHMRHist *histos;

    std::vector<std::vector<int>> photonId; // IDs of the photon pairs for pi0 candidates
    std::vector<std::vector<int>> sidebandId; // IDs of the photon pairs for sideband candidates
};

#endif
