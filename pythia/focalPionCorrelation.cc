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

void DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins]);
void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins]);
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

    TF1 *fPhotonEfficiency = new TF1("fPhotonEfficiency", "TMath::Exp(-3.20093/x)"); // Parameters from fit to efficiency (PhotonEfficiency.C)
    TRandom3 *rand = new TRandom3();

    // Basic histograms
    TH1D *hCounter = new TH1D("hCounter", "hCounter", 10, 0, 10);
    TH1D *hRealTriggCounter = new TH1D("hRealTriggCounter", "hRealtriggCounter", 10, 0, 10);
    
    for (int i = 0; i <= nIncPtBin; i++) logBinsX[i] = limMin*exp(i*logBW);
    
    TH1D *hPionPt = new TH1D("hPionPt", "hPionPt", nIncPtBin, logBinsX); hPionPt->Sumw2();
    TH1D *hPionPtFor = new TH1D("hPionPtFor", "hPionPtFor", nIncPtBin, logBinsX); hPionPtFor->Sumw2();
    TH1D *hPionPtMid = new TH1D("hPionPtMid", "hPionPtMid", nIncPtBin, logBinsX); hPionPtMid->Sumw2();
    
    TH1D *hChargedHadronPt = new TH1D("hChargedHadronPt", "hChargedHadronPt", nIncPtBin, logBinsX); hChargedHadronPt->Sumw2();
    TH1D *hChargedHadronPtFor = new TH1D("hChargedHadronPtFor", "hChargedHadronPtFor", nIncPtBin, logBinsX); hChargedHadronPtFor->Sumw2();
    TH1D *hChargedHadronPtMid = new TH1D("hChargedHadronPtMid", "hChargedHadronPtMid", nIncPtBin, logBinsX); hChargedHadronPtMid->Sumw2();

    TH1D *hPhotonPt = new TH1D("hPhotonPt", "hPhotonPt", nIncPtBin, logBinsX); hPhotonPt->Sumw2();
    TH1D *hPhotonPtFor = new TH1D("hPhotonPtFor", "hPhotonPtFor", nIncPtBin, logBinsX); hPhotonPtFor->Sumw2();
    TH1D *hPhotonPtMid = new TH1D("hPhotonPtMid", "hPhotonPtMid", nIncPtBin, logBinsX); hPhotonPtMid->Sumw2();
    TH1D *hPhotonEnergyReal = new TH1D("hPhotonEnergyReal", "hPhotonEnergyReal", nPhotonEnergyBin, limPhotonEnergyMin, limPhotonEnergyMax); hPhotonEnergyReal->Sumw2();
    TH1D *hPhotonEnergy = new TH1D("hPhotonEnergy", "hPhotonEnergy", nPhotonEnergyBin, limPhotonEnergyMin, limPhotonEnergyMax); hPhotonEnergy->Sumw2();
    
    TH1D *hPionEta = new TH1D("hPionEta", "hPionEta", nIncEtaBin, -incEtaRange/2., incEtaRange/2.); hPionEta->Sumw2();
    TH1D *hChargedHadronEta = new TH1D("hChargedHadronEta", "hChargedHadronEta", nIncEtaBin, -incEtaRange/2., incEtaRange/2.); hChargedHadronEta->Sumw2();

    // Correlation and mass histograms
    TDirectory *dirMasses = fOut->mkdir("Masses");
    TDirectory *dirCorrMid = fOut->mkdir("CorrMid");
    TDirectory *dirCorrFor = fOut->mkdir("CorrFor");
    TDirectory *dirCorrChargedMid = fOut->mkdir("CorrChargedMid");
    TDirectory *dirCorrChargedFor = fOut->mkdir("CorrChargedFor");
    TDirectory *dirCorrMassMass = fOut->mkdir("CorrMassMass");
    TDirectory *dirCorrMassSide = fOut->mkdir("CorrMassSide");
    TDirectory *dirCorrSideMass = fOut->mkdir("CorrSideMass");
    TDirectory *dirCorrSideSide = fOut->mkdir("CorrSideSide");

    TH1D *hPi0MassTrigg[nTriggBins];
    TH1D *hPi0MassAssocPeak[nTriggBins][nAssocBins];
    TH1D *hPi0MassAssocSide[nTriggBins][nAssocBins];
 
    TH2D *hCorrMid[nTriggBins][nAssocBins];
    TH2D *hCorrFor[nTriggBins][nAssocBins];
    TH2D *hCorrChargedMid[nTriggBins][nAssocBins];
    TH2D *hCorrChargedFor[nTriggBins][nAssocBins];
    
    TH2D *hCorrMassMass[nTriggBins][nAssocBins];
    TH2D *hCorrMassSide[nTriggBins][nAssocBins];
    TH2D *hCorrSideMass[nTriggBins][nAssocBins];
    TH2D *hCorrSideSide[nTriggBins][nAssocBins];

    for (int i = 0; i < nTriggBins; i++) {
        double tlow = triggPt[i];
        double tupp = triggPt[i+1];
        dirMasses->cd();
        hPi0MassTrigg[i] = new TH1D(Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), 360, 0.0, 720.0);
        hPi0MassTrigg[i]->Sumw2();

        for (int j = 0; j < nAssocBins; j++) {
            double alow = assocPt[j];
            double aupp = assocPt[j+1];

            if (tlow < aupp) continue;

            dirMasses->cd();
            hPi0MassAssocPeak[i][j] = new TH1D(Form("hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                        Form("hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 360, 0.0, 720.0);
            hPi0MassAssocPeak[i][j]->Sumw2();
            hPi0MassAssocSide[i][j] = new TH1D(Form("hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                        Form("hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 360, 0.0, 720.0);
            hPi0MassAssocSide[i][j]->Sumw2();

            dirCorrMid->cd();
            hCorrMid[i][j] = new TH2D(Form("hCorrMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -2.*etaTrackerRange, 2.*etaTrackerRange);
            hCorrMid[i][j]->Sumw2();
            
            dirCorrFor->cd();
            hCorrFor[i][j] = new TH2D(Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrFor[i][j]->Sumw2();
            
            dirCorrChargedMid->cd();
            hCorrChargedMid[i][j] = new TH2D(Form("hCorrChargedMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrChargedMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -2.*etaTrackerRange, 2.*etaTrackerRange);
            hCorrChargedMid[i][j]->Sumw2();
            
            dirCorrChargedFor->cd();
            hCorrChargedFor[i][j] = new TH2D(Form("hCorrChargedFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrChargedFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrChargedFor[i][j]->Sumw2();
            
            dirCorrMassMass->cd();
            hCorrMassMass[i][j] = new TH2D(Form("hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrMassMass[i][j]->Sumw2();

            dirCorrMassSide->cd();
            hCorrMassSide[i][j] = new TH2D(Form("hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrMassSide[i][j]->Sumw2();

            dirCorrSideMass->cd();
            hCorrSideMass[i][j] = new TH2D(Form("hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrSideMass[i][j]->Sumw2();

            dirCorrSideSide->cd();
            hCorrSideSide[i][j] = new TH2D(Form("hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrSideSide[i][j]->Sumw2();
            
        }
    }

    fOut->cd();
    
    // Particle lists
    TClonesArray *arrPi0Mid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPi0For = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrChargedMid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrChargedFor = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPhotonMid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPhotonFor = new TClonesArray("TLorentzVector", 1500);

    TClonesArray *arrPi0Peak = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPi0Side = new TClonesArray("TLorentzVector", 1500);
 
    // 
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        hCounter->Fill(0.5); // Number of events

        std::cout << "event " << iEvent << std::endl;

        if ( !pythia.next() ) continue;

        int nPi0Mid = 0;
        int nPi0For = 0;
        int nChargedMid = 0;
        int nChargedFor = 0;
        int nPhotonMid = 0;
        int nPhotonFor = 0;

        arrPi0Mid->Clear();
        arrPi0For->Clear();
        arrChargedMid->Clear();
        arrChargedFor->Clear();
        arrPhotonMid->Clear();
        arrPhotonFor->Clear();
    
        arrPi0Peak->Clear(); 
        arrPi0Side->Clear();

        // Collect particles of interest
        for (int iPart = 0; iPart < pythia.event.size(); iPart++) {

            TLorentzVector lv(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), pythia.event[iPart].e());

            if (pythia.event[iPart].id()==22) { // photons
                double px = pythia.event[iPart].px();
                double py = pythia.event[iPart].py();
                double pz = pythia.event[iPart].pz();
                double e = TMath::Sqrt(px*px + py*py + pz*pz);
                double eSmear = PhotonEnergySmearing(rand, px, py, pz);
     
                hPhotonEnergyReal->Fill(eSmear);

                if (!IsPhotonRemoved(eSmear, rand, fPhotonEfficiency)) {
                    hPhotonEnergy->Fill(eSmear);
                    TLorentzVector lvSmeared = TLorentzVector(eSmear*px/e, eSmear*py/e, eSmear*pz/e, eSmear);
               
                    hPhotonPt->Fill(TMath::Sqrt(px*px + py*py)); 
                    double eta = pythia.event[iPart].eta();
                    if (IsTrackerAcceptance(eta)) {
                        new((*arrPhotonMid)[nPhotonMid++]) TLorentzVector(lvSmeared);
                        hPhotonPtMid->Fill(TMath::Sqrt(px*px + py*py)); 
                    }
                    if (IsFocalAcceptance(eta)) {
                        new((*arrPhotonFor)[nPhotonFor++]) TLorentzVector(lvSmeared);
                        hPhotonPtFor->Fill(TMath::Sqrt(px*px + py*py)); 
                    }
                }
            }

            if (pythia.event[iPart].id()==111) { // Pions
                double pt = pythia.event[iPart].pT();
                double eta = pythia.event[iPart].eta();

                hPionPt->Fill(pt);
                hPionEta->Fill(eta);
                if (IsTrackerAcceptance(eta)) {
                    hPionPtMid->Fill(pt);
                    new((*arrPi0Mid)[nPi0Mid++]) TLorentzVector(lv);
                }
                if (IsFocalAcceptance(eta)) {
                    hPionPtFor->Fill(pt);
                    new((*arrPi0For)[nPi0For++]) TLorentzVector(lv);
                }
            }

            if (pythia.event[iPart].isFinal() && pythia.event[iPart].isHadron() && pythia.event[iPart].isCharged()) { // charged hadrons
                double pt = pythia.event[iPart].pT();
                double eta = pythia.event[iPart].eta();

                hChargedHadronPt->Fill(pt);
                hChargedHadronEta->Fill(eta);
                if (IsTrackerAcceptance(eta)) {
                    hChargedHadronPtMid->Fill(pt);
                    new((*arrChargedMid)[nChargedMid++]) TLorentzVector(lv);
               }
                if (IsFocalAcceptance(eta)) {
                    hChargedHadronPtFor->Fill(pt);
                    new((*arrChargedFor)[nChargedFor++]) TLorentzVector(lv);
                }
            }
        }

        if (nPi0Mid>0) hCounter->Fill(1.5); // number of events with pion0 in mid rapidity
        if (nPi0For>0) hCounter->Fill(2.5); // number of events with pion0 in forward rapidity 
        
        ReconstructPions(arrPhotonFor, arrPi0Peak, 1);
        ReconstructPions(arrPhotonFor, arrPi0Side, 0);

        std::vector<int> listTriggReal, listAssocReal, listTriggPeak, listTriggSide, listAssocPeak, listAssocSide;
        int binsWithTriggReal[nTriggBins+1] = {0}, binsWithTriggPeak[nTriggBins+1] = {0}, binsWithTriggSide[nTriggBins+1] = {0};
        GetTriggAssocLists(arrPi0For, listTriggReal, listAssocReal, binsWithTriggReal, bUseLeading); 
        GetTriggAssocLists(arrPi0Peak, listTriggPeak, listAssocPeak, binsWithTriggPeak, bUseLeading); 
        GetTriggAssocLists(arrPi0Side, listTriggSide, listAssocSide, binsWithTriggSide, bUseLeading); 
    
        FillPionMasses(arrPhotonFor, hPi0MassTrigg, hPi0MassAssocPeak, hPi0MassAssocSide, binsWithTriggPeak, binsWithTriggSide);
        FillRealTriggers(hRealTriggCounter, arrPi0For, listTriggReal);
       
        DoCorrelations(arrPi0For, listTriggReal, listAssocReal, hCorrFor);
        DoCorrelations(arrPi0Peak, listTriggPeak, listAssocPeak, hCorrMassMass);
        DoCorrelations(arrPi0Side, listTriggSide, listAssocSide, hCorrSideSide);
        if (bUseLeading) {
            int isPeakTriggLarger = GetLargerTrigg(arrPi0Peak, listTriggPeak, arrPi0Side, listTriggSide);
            if (isPeakTriggLarger) {
                DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, hCorrMassSide);
            } else {
                DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, hCorrSideMass);
            }
        } else {
            DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, hCorrMassSide);
            DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, hCorrSideMass);
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

void DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins])
{
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

        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];
            if (iTrigg==iAssoc) continue; // autocorrelations
            TLorentzVector *lvAssoc = (TLorentzVector*)arrPi0->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, nAssocBins, ptAssoc);

            if (triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;
            
            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;
            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta);
        }
    }

}

void DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins])
{
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

        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];

            TLorentzVector *lvAssoc = (TLorentzVector*)arrPi0Assoc->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, nAssocBins, ptAssoc);

            if (triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;
            
            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;
            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta);
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
