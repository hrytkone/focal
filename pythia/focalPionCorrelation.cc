#include <vector>
#include <cstring>

#include "Pythia8/Pythia.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TStopwatch.h"

#include "AliJHMRHist.h"

using namespace Pythia8;

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

//-------------------------
//        Functions       |
//-------------------------
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
void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins], bool bUseWeightTrigg, bool bUseWeightAssoc, TF1 *fPhotonAcceptanceEfficiency);
void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, bool bMass);
void GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bUseLeading);
void FillRealTriggers(TH1D *hRealTriggCounter, TClonesArray *arrRealPi0, std::vector<int>& listTrigg);
void FillPionMasses(TClonesArray *arrPhoton, TH1D *hMassTrigg[nTriggBins], TH1D *hMassAssocPeak[nTriggBins][nAssocBins], TH1D *hMassAssocSide[nTriggBins][nAssocBins], int binsWithTriggPeak[nTriggBins], int binsWithTriggSide[nTriggBins]);

int main(int argc, char *argv[]) {

    if (argc==1) {
        cout << "Usage : ./focalPionCorrelation output.root bUseLeading pythiaSettings.cmnd seed" << endl;
        return 0;
    }

    TStopwatch timer;
    timer.Start();

    TString outFileName = argc > 1 ? argv[1] : "output.root";
    int bUseLeading = argc > 2 ? atol(argv[2]) : 0;
    TString pythiaSettings = argc > 3 ? argv[3] : "PythiaHard.cmnd";
    int seed = argc > 4 ? atol(argv[4]) : 0;

    TFile *fOut = new TFile(outFileName, "RECREATE");

    Pythia pythia;

    // Initialise pythia
    pythia.readFile(pythiaSettings.Data());
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.init();

    int nEvents = pythia.mode("Main:numberOfEvents");

    JHMRHist *fHistos = new JHMRHist();
    JHMRPythiaCatalyst *fCatalyst = new JHMRPythiaCatalyst(histos);
    JHMRCorr *fCorr = new JHMRCorr(histos);

    TF1 *fPhotonEfficiency = new TF1("fPhotonEfficiency", "TMath::Exp(-3.20093/x)"); // Parameters from fit to efficiency (PhotonEfficiency.C)
    TF1 *fPhotonAcceptanceEfficiency = new TF1("fPhotonAcceptanceEfficiency", "TMath::Exp(-0.117082/(x + 0.0832931))"); // Parameters from fit (CheckMissingPionsRatio.C)
    TRandom3 *rand = new TRandom3();

    fOut->cd();

    // 
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        //std::cout << "event " << iEvent << std::endl;

        if ( !pythia.next() ) continue;
        
        fCatalyst->ReconstructPions(arrPhotonFor, arrPi0Peak, 1);
        fCatalyst->ReconstructPions(arrPhotonFor, arrPi0Side, 0);

        std::vector<int> listTriggReal, listAssocReal, listTriggPeak, listTriggSide, listAssocPeak, listAssocSide;
        int binsWithTriggReal[nTriggBins+1] = {0}, binsWithTriggPeak[nTriggBins+1] = {0}, binsWithTriggSide[nTriggBins+1] = {0};
        fCorr->GetTriggAssocLists(arrPi0For, listTriggReal, listAssocReal, binsWithTriggReal, bUseLeading); 
        fCorr->GetTriggAssocLists(arrPi0Peak, listTriggPeak, listAssocPeak, binsWithTriggPeak, bUseLeading); 
        fCorr->GetTriggAssocLists(arrPi0Side, listTriggSide, listAssocSide, binsWithTriggSide, bUseLeading); 
    
        FillPionMasses(arrPhotonFor, hPi0MassTrigg, hPi0MassAssocPeak, hPi0MassAssocSide, binsWithTriggPeak, binsWithTriggSide);
        FillRealTriggers(hRealTriggCounter, arrPi0For, listTriggReal);
       
        fCorr->DoCorrelations(arrPi0For, listTriggReal, listAssocReal, hCorrFor, 0, fPhotonAcceptanceEfficiency);
        fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, listAssocPeak, hCorrMassMass, 1, fPhotonAcceptanceEfficiency);
        fCorr->DoCorrelations(arrPi0Side, listTriggSide, listAssocSide, hCorrSideSide, 0, fPhotonAcceptanceEfficiency);
        if (bUseLeading) {
            int isPeakTriggLarger = fCorr->GetLargerTrigg(arrPi0Peak, listTriggPeak, arrPi0Side, listTriggSide);
            if (isPeakTriggLarger) {
                fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, hCorrMassSide, 1, 0, fPhotonAcceptanceEfficiency);
            } else {
                fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, hCorrSideMass, 0, 1, fPhotonAcceptanceEfficiency);
            }
        } else {
            fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, hCorrMassSide, 1, 0, fPhotonAcceptanceEfficiency);
            fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, hCorrSideMass, 0, 1, fPhotonAcceptanceEfficiency);
        }
    }

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    pythia.stat();

    timer.Print();

    return 0;
}

//-------------------------
//        Functions       |
//-------------------------

bool IsTrackerAcceptance(double eta, double etaRange)
{
    return TMath::Abs(eta) < etaRange/2. ? true : false;
}

bool IsFocalAcceptance(double eta, double etaMin, double etaMax)
{
    return (eta > etaMin && eta < etaMax) ? true : false;
}

bool IsMassWindow(double mass)
{
    return (mass > 110. && mass < 160.) ? true : false;
}

bool IsSideband(double mass)
{
    return ((mass > 40. && mass < 80.) || (mass > 210. && mass < 280.)) ? true : false;
}

// Return 1 if trigger from peak has larger pT than trigger from sideband, 0 if sideband
// trigger pT is larger
int GetLargerTrigg(TClonesArray *arrPi0Peak, std::vector<int> listTriggPeak, TClonesArray *arrPi0Side, std::vector<int> listTriggSide)
{
    if (!listTriggPeak.size()) return 0;
    if (!listTriggSide.size()) return 1;
    TLorentzVector *lvTriggPeak = (TLorentzVector*)arrPi0Peak->At(listTriggPeak[0]);
    TLorentzVector *lvTriggSide = (TLorentzVector*)arrPi0Side->At(listTriggSide[0]);
    double pTpeak = lvTriggPeak->Pt();
    double pTside = lvTriggSide->Pt();
    if (pTpeak > pTside) return 1;
    return 0;
}

void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, bool bMass)
{
    int nCandidate = 0;
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector lvSum = GetPhotonSumVector(arrPhoton, i, j);
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            bool bIsInWindow = bMass ? IsMassWindow(mass) : IsSideband(mass); 
            if (bIsInWindow) new ((*arrPi0Candidates)[nCandidate++]) TLorentzVector(lvSum);
        }
    }
}

void FillRealTriggers(TH1D *hRealTriggCounter, TClonesArray *arrRealPi0, std::vector<int>& listTrigg)
{
    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;
    
    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrRealPi0->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        int iTriggBin = GetBin(triggPt, nTriggBins, ptTrigg);
        hRealTriggCounter->Fill(iTriggBin+0.5);
    }
}

void GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bUseLeading)
{
    int iLeadingTrigg = -1;
    if (bUseLeading) {
        iLeadingTrigg = GetLeadingTriggerIndex(arrPi0Candidates);
        if (iLeadingTrigg > -1) listTrigg.push_back(iLeadingTrigg);
    }

    int nPi0Candidates = arrPi0Candidates->GetEntriesFast(); 
    for (int i = 0; i < nPi0Candidates; i++) {
        TLorentzVector *lvPion = (TLorentzVector*)arrPi0Candidates->At(i);
        double pT = lvPion->Pt();
        int iTrigg = GetBin(triggPt, nTriggBins, pT);
        int iAssoc = GetBin(assocPt, nAssocBins, pT);
        if (bUseLeading) {
            if (i != iLeadingTrigg && iAssoc >= 0) listAssoc.push_back(i);
        } else {
            if (iTrigg >= 0) listTrigg.push_back(i);
            if (iAssoc >= 0) listAssoc.push_back(i);
        }
        if (iTrigg > 0) binsWithTrigg[iTrigg]++;
    }
}

void DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins], bool bUseWeight, TF1 *fPhotonAcceptanceEfficiency)
{
    double wTrigg = 1.0;
    double wAssoc = 1.0;

    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;
    
    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrPi0->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();
        
        int iTriggBin = GetBin(triggPt, nTriggBins, ptTrigg);
               
        int nAssoc = listAssoc.size();       
        if (nAssoc < 1) continue;

        if (bUseWeight) wTrigg = 1./fPhotonAcceptanceEfficiency->Eval(ptTrigg);
        
        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];
            if (iTrigg==iAssoc) continue; // autocorrelations
            TLorentzVector *lvAssoc = (TLorentzVector*)arrPi0->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, nAssocBins, ptAssoc);

            if (triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;
            
            if (bUseWeight) wAssoc = 1./fPhotonAcceptanceEfficiency->Eval(ptAssoc); 

            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;
            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
        }
    }

}

void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins], bool bUseWeightTrigg, bool bUseWeightAssoc, TF1 *fPhotonAcceptanceEfficiency)
{
    double wTrigg = 1.0;
    double wAssoc = 1.0;
    
    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrPi0Trigg->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();
        int iTriggBin = GetBin(triggPt, nTriggBins, ptTrigg);
        
        int nAssoc = listAssoc.size();       
        if (nAssoc < 1) continue;
        
        if (bUseWeightTrigg) wTrigg = 1./fPhotonAcceptanceEfficiency->Eval(ptTrigg);

        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];

            TLorentzVector *lvAssoc = (TLorentzVector*)arrPi0Assoc->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, nAssocBins, ptAssoc);

            if (triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;
            
            if (bUseWeightAssoc) wAssoc = 1./fPhotonAcceptanceEfficiency->Eval(ptAssoc);
            
            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;
            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
        }
    }
}

int GetBin(double arr[], int nArr, double val)
{   
    for (int i=0; i<nArr; i++) {
        if (arr[i]<=val && val<arr[i+1]) return i;
    }
    return -1;
}

double GetDeltaPhi(double phiTrigg, double phiAssoc)
{
    double dphi = phiTrigg - phiAssoc;
    if (dphi>3.0*TMath::Pi()/2.0) return dphi - 2.0*TMath::Pi();
    if (dphi<-TMath::Pi()/2.0) return dphi + 2.0*TMath::Pi();
    return dphi;
}

// For FoCal:
//      sigma/E = a/sqrt(E) + b = 0.27/sqrt(E) + 0.01
double PhotonEnergySmearing(TRandom3 *rand, double px, double py, double pz)
{
    const double eCut = 0.01;
    double e = sqrt(px*px + py*py + pz*pz);
    if (e < eCut) return e;

    double a = 0.27;
    double b = 0.01;
    double sigma = sqrt(a*a*e + b*b*e*e);
    double eSmear = -1.0;
    while (eSmear < 0.0) eSmear = rand->Gaus(e, sigma);
    return eSmear;
}

// To count in single photon efficiency
bool IsPhotonRemoved(double ePhoton, TRandom3 *rand, TF1 *fPhotonEff)
{
    double photonEff = fPhotonEff->Eval(ePhoton);
    if (rand->Uniform() > photonEff) return true;
    return false;
}

TLorentzVector GetPhotonSumVector(TClonesArray *arrPhoton, int iPhoton1, int iPhoton2)
{
    TLorentzVector *lv1 = (TLorentzVector*)arrPhoton->At(iPhoton1);
    TLorentzVector *lv2 = (TLorentzVector*)arrPhoton->At(iPhoton2);
    TLorentzVector lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}

int GetLeadingTriggerIndex(TClonesArray *arrPi0)
{
    int itrigg = -1;
    double pTmax = 0;
    int nPi0 = arrPi0->GetEntriesFast();
    for (int i = 0; i < nPi0; i++) {
        TLorentzVector *lvPi0 = (TLorentzVector*)arrPi0->At(i);
        double pT = lvPi0->Pt();
        if (pT > triggPt[0] && pTmax < pT) {
            pTmax = pT;
            itrigg = i;
        }
    }
    return itrigg;
}


void FillPionMasses(TClonesArray *arrPhoton, TH1D *hMassTrigg[nTriggBins], TH1D *hMassAssocPeak[nTriggBins][nAssocBins], TH1D *hMassAssocSide[nTriggBins][nAssocBins], int binsWithTriggPeak[nTriggBins], int binsWithTriggSide[nTriggBins])
{
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector lvSum = GetPhotonSumVector(arrPhoton, i, j);
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            int iTriggBin = GetBin(triggPt, nTriggBins, pT);
            int iAssocBin = GetBin(assocPt, nAssocBins, pT);
            if (iTriggBin >= 0) hMassTrigg[iTriggBin]->Fill(mass);
            for (int it = 0; it < nTriggBins; it++) {
                if (triggPt[it] < assocPt[iAssocBin+1]) continue;
                if (binsWithTriggPeak[it] > 0 && iAssocBin >= 0) hMassAssocPeak[it][iAssocBin]->Fill(mass);
                if (binsWithTriggSide[it] > 0 && iAssocBin >= 0) hMassAssocSide[it][iAssocBin]->Fill(mass);
            }
        }
    }
}
