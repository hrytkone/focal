#ifndef ALIJHMRDATAMANAGER_H
#define ALIJHMRDATAMANAGER_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <Pythia8/Pythia.h>

#include <TString.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include <vector>

using namespace std;
using namespace Pythia8;

class AliJHMRDataManager {
  
public:
    
    AliJHMRDataManager():

    bool IsTrackerAcceptance(double eta, double etaRange=etaTrackerRange);
    bool IsFocalAcceptance(double eta, double etaMin=etaFocalMin, double etaMax=etaFocalMax);
    int GetBin(double arr[], int nArr, double val);
    double GetDeltaPhi(double phiTrigg, double phiAssoc);
    double PhotonEnergySmearing(TRandom3 *rand, double px, double py, double pz);
    bool IsPhotonRemoved(double ePhoton, TRandom3 *rand, TF1 *fPhotonEff);
TLorentzVector GetPhotonSumVector(TClonesArray *arrayPhoton, int iPhoton1, int iPhoton2);
    int GetLeadingTriggerIndex(TClonesArray *arrPi0);
    int GetLargerTrigg(TClonesArray *arrPi0Peak, std::vector<int> listTriggPeak, TClonesArray *arrPi0Side, std::vector<int> listTriggSide);

    void DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins], bool bUseWeight, TF1 *fPhotonAcceptanceEfficiency);
    void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins], bool bUseWeightTrigg, bool bUseWeightAssoc,  TF1 *fPhotonAcceptanceEfficiency);
    void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, bool bMass);
    void GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bUseLeading);

protected:

    //-------------------------
    //         pT bins        |
    //-------------------------
    const int nTriggBins = 4;
    double  triggPt[nTriggBins+1] = {1.0, 2.0, 4.0, 8.0, 20.0};
    //double  triggPt[nTriggBins+1] = {3.0, 1000.0};
    //double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    const int nAssocBins = 4;
    double  assocPt[nAssocBins+1] = {0.5, 1.0, 2.0, 3.0, 4.0};
    //double  assocPt[nAssocBins+1] = {3.0, 1000.0};
    //double  assocPt[nAssocBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};
  
};

#endif

