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

using namespace Pythia8;

//-------------------------
//         pT bins        |
//-------------------------
const int nTriggBins = 8;
double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

const int nAssocBins = 7;
double  assocPt[nAssocBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};

const int nPi0PtBins = 15;
double pi0Pt[nPi0PtBins+1] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 50.0};

//-------------------------
//        Functions       |
//-------------------------
void DoCorrelations(TClonesArray *arrPions, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins]);
int GetBin(double arr[], int nArr, double val);
void GetTriggAssocLists(TClonesArray *arrPions, std::vector<int>& listTrigg, std::vector<int>& listAssoc);
double GetDeltaPhi(double phiTrigg, double phiAssoc);

double PhotonEnergySmearing(TRandom3 *rand, double px, double py, double pz);
void ReconstructPions(TClonesArray *arrayPhoton, TClonesArray *arrayPi0Rec, TClonesArray *arrayPi0Side);
void FillPionMasses(TClonesArray *arrayPhoton, TH1D *hMasses[nPi0PtBins]);

int main(int argc, char *argv[]) {

    if (argc==1) {
        cout << "Usage : ./focalPionCorrelation output.root pythiaSettings.cmnd seed" << endl;
        return 0;
    }

    TString outFileName = argc > 1 ? argv[1] : "output.root";
    TString pythiaSettings = argc > 2 ? argv[2] : "PythiaHard.cmnd";
    int seed = argc > 3 ? atol(argv[3]) : 0;

    TFile *fOut = new TFile(outFileName, "RECREATE");

    Pythia pythia;

    // Initialise pythia
    pythia.readFile(pythiaSettings.Data());
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.init();

    int nEvents = pythia.mode("Main:numberOfEvents");

    TRandom3 *rand = new TRandom3();

    // Basic histograms
    TH1D *hCounter = new TH1D("hCounter", "hCounter", 10, 0, 10);
    
    TH1D *hPionPt = new TH1D("hPionPt", "hPionPt", 401, 0.0, 15.0); hPionPt->Sumw2();
    TH1D *hPionPtFor = new TH1D("hPionPtFor", "hPionPtFor", 401, 0.0, 15.0); hPionPtFor->Sumw2();
    TH1D *hPionPtMid = new TH1D("hPionPtMid", "hPionPtMid", 401, 0.0, 15.0); hPionPtMid->Sumw2();
    
    TH1D *hChargedHadronPt = new TH1D("hChargedHadronPt", "hChargedHadronPt", 401, 0.0, 15.0); hChargedHadronPt->Sumw2();
    TH1D *hChargedHadronPtFor = new TH1D("hChargedHadronPtFor", "hChargedHadronPtFor", 401, 0.0, 15.0); hChargedHadronPtFor->Sumw2();
    TH1D *hChargedHadronPtMid = new TH1D("hChargedHadronPtMid", "hChargedHadronPtMid", 401, 0.0, 15.0); hChargedHadronPtMid->Sumw2();

    TH1D *hPionEta = new TH1D("hPionEta", "hPionEta", 400, -15.0, 15.0); hPionEta->Sumw2();
    TH1D *hChargedHadronEta = new TH1D("hChargedHadronEta", "hChargedHadronEta", 401, -15.0, 15.0); hChargedHadronEta->Sumw2();

    // Correlation histograms
    const int NphiBins = int( 6.3/0.05 ) + 1;
    const double etaRange = 2.0;
    const double etaBinWidth = 0.025e0;
    const int NetaBins = int( 2.2*etaRange / etaBinWidth ) + 1;

    TH2D *hCorrMid[nTriggBins][nAssocBins];
    TH2D *hCorrFor[nTriggBins][nAssocBins];
    TH2D *hCorrChargedMid[nTriggBins][nAssocBins];
    TH2D *hCorrChargedFor[nTriggBins][nAssocBins];
    TH2D *hCorrPion0RecMid[nTriggBins][nAssocBins];
    TH2D *hCorrPion0RecFor[nTriggBins][nAssocBins];
    TH2D *hCorrPion0SideMid[nTriggBins][nAssocBins];
    TH2D *hCorrPion0SideFor[nTriggBins][nAssocBins];
    for (int i = 0; i < nTriggBins; i++) {
        for (int j = 0; j < nAssocBins; j++) {
            hCorrMid[i][j] = new TH2D(Form("hCorrMid%d:%d", i, j), Form("hCorrMid%d:%d", i, j), NphiBins, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, NetaBins, -3.0, 3.0);
            hCorrMid[i][j]->Sumw2();
            hCorrFor[i][j] = new TH2D(Form("hCorrFor%d:%d", i, j), Form("hCorrFor%d:%d", i, j), NphiBins, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, NetaBins, -3.0, 3.0);
            hCorrFor[i][j]->Sumw2();
            hCorrChargedMid[i][j] = new TH2D(Form("hCorrChargedMid%d:%d", i, j), Form("hCorrChargedMid%d:%d", i, j), NphiBins, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, NetaBins, -3.0, 3.0);
            hCorrChargedMid[i][j]->Sumw2();
            hCorrChargedFor[i][j] = new TH2D(Form("hCorrChargedFor%d:%d", i, j), Form("hCorrChargedFor%d:%d", i, j), NphiBins, -TMath::Pi()/4.0, 3.0*TMath::Pi()/4.0, NetaBins, -3.0, 3.0);
            hCorrChargedFor[i][j]->Sumw2();
            hCorrPion0RecMid[i][j] = new TH2D(Form("hCorrPion0RecMid%d:%d", i, j), Form("hCorrPion0RecMid%d:%d", i, j), NphiBins, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, NetaBins, -3.0, 3.0);
            hCorrPion0RecMid[i][j]->Sumw2();
            hCorrPion0RecFor[i][j] = new TH2D(Form("hCorrPion0RecFor%d:%d", i, j), Form("hCorrPion0RecFor%d:%d", i, j), NphiBins, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, NetaBins, -3.0, 3.0);
            hCorrPion0RecFor[i][j]->Sumw2();
            hCorrPion0SideMid[i][j] = new TH2D(Form("hCorrPion0SideMid%d:%d", i, j), Form("hCorrPion0SideMid%d:%d", i, j), NphiBins, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, NetaBins, -3.0, 3.0);
            hCorrPion0SideMid[i][j]->Sumw2();
            hCorrPion0SideFor[i][j] = new TH2D(Form("hCorrPion0SideFor%d:%d", i, j), Form("hCorrPion0SideFor%d:%d", i, j), NphiBins, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, NetaBins, -3.0, 3.0);
            hCorrPion0SideFor[i][j]->Sumw2();
        }
    }

    // Particle lists
    TClonesArray *arrPion0Mid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0For = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrChargedMid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrChargedFor = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPhotonMid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPhotonFor = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0RecMid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0RecFor = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0SideMid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0SideFor = new TClonesArray("TLorentzVector", 1500);

    //
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        hCounter->Fill(0.5); // Number of events

        if ( !pythia.next() ) continue;

        int nPion0Mid = 0;
        int nPion0For = 0;
        int nChargedMid = 0;
        int nChargedFor = 0;
        int nPhotonMid = 0;
        int nPhotonFor = 0;

        arrPion0Mid->Clear();
        arrPion0For->Clear();
        arrChargedMid->Clear();
        arrChargedFor->Clear();
        arrPhotonMid->Clear();
        arrPhotonFor->Clear();
        arrPion0RecFor->Clear(); 
        arrPion0RecMid->Clear(); 
        arrPion0SideFor->Clear(); 
        arrPion0SideMid->Clear(); 

        // Trigg & associated lists, store indices
        std::vector<int> listTriggMid, listTriggFor, listTriggChargedMid, listTriggChargedFor, listTriggPion0RecMid, listTriggPion0RecFor, listTriggPion0SideMid, listTriggPion0SideFor;
        std::vector<int> listAssocMid, listAssocFor, listAssocChargedMid, listAssocChargedFor, listAssocPion0RecMid, listAssocPion0RecFor, listAssocPion0SideMid, listAssocPion0SideFor;

        // Collect particles of interest
        for (int iPart = 0; iPart < pythia.event.size(); iPart++) {

            TLorentzVector lv(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), pythia.event[iPart].e());

            if (pythia.event[iPart].id()==22) { // photons
                double px = pythia.event[iPart].px();
                double py = pythia.event[iPart].py();
                double pz = pythia.event[iPart].pz();
                double e = TMath::Sqrt(px*px + py*py + pz*pz);
                double eSmear = PhotonEnergySmearing(rand, px, py, pz);
                TLorentzVector lv = TLorentzVector(eSmear*px/e, eSmear*py/e, eSmear*pz/e, eSmear);
                
                double eta = pythia.event[iPart].eta();
                if (TMath::Abs(eta)<1.0) {
                    new((*arrPhotonMid)[nPhotonMid++]) TLorentzVector(lv);
                }
                if (eta>3.2 && eta<5.8) {
                    new((*arrPhotonFor)[nPhotonFor++]) TLorentzVector(lv);
                }

            }

            if (pythia.event[iPart].id()==111) { // Pions
                double pt = pythia.event[iPart].pT();
                double eta = pythia.event[iPart].eta();

                hPionPt->Fill(pt);
                hPionEta->Fill(eta);
                if (TMath::Abs(eta)<1.0) {
                    hPionPtMid->Fill(pt);
                    new((*arrPion0Mid)[nPion0Mid++]) TLorentzVector(lv);
                }
                if (eta>3.2 && eta<5.8) {
                    hPionPtFor->Fill(pt);
                    new((*arrPion0For)[nPion0For++]) TLorentzVector(lv);
                }
            }

            if (pythia.event[iPart].isFinal() && pythia.event[iPart].isHadron() && (pythia.event[iPart].charge() > 0)) { // charged hadrons
                double pt = pythia.event[iPart].pT();
                double eta = pythia.event[iPart].eta();

                hChargedHadronPt->Fill(pt);
                hChargedHadronEta->Fill(eta);
                if (TMath::Abs(eta)<1.0) {
                    hChargedHadronPtMid->Fill(pt);
                    new((*arrChargedMid)[nChargedMid++]) TLorentzVector(lv);
                }
                if (eta>3.2 && eta<5.8) {
                    hChargedHadronPtFor->Fill(pt);
                    new((*arrChargedFor)[nChargedFor++]) TLorentzVector(lv);
                }
            }
        }

        ReconstructPions(arrPhotonMid, arrPion0RecMid, arrPion0SideMid);
        ReconstructPions(arrPhotonFor, arrPion0RecFor, arrPion0SideFor);

        if (nPion0Mid>0) hCounter->Fill(1.5); // number of events with pion0 in mid rapidity
        if (nPion0For>0) hCounter->Fill(2.5); // number of events with pion0 in forward rapidity 
        if (nChargedMid>0) hCounter->Fill(3.5); // number of events with charged hadron in mid rapidity
        if (nChargedFor>0) hCounter->Fill(4.5); // number of events with charged hadron in forward rapidity

        GetTriggAssocLists(arrPion0Mid, listTriggMid, listAssocMid);
        GetTriggAssocLists(arrPion0For, listTriggFor, listAssocFor);
        GetTriggAssocLists(arrChargedFor, listTriggChargedFor, listAssocChargedFor);
        GetTriggAssocLists(arrChargedMid, listTriggChargedMid, listAssocChargedMid);
        GetTriggAssocLists(arrPion0RecMid, listTriggPion0RecMid, listAssocPion0RecMid);
        GetTriggAssocLists(arrPion0RecFor, listTriggPion0RecFor, listAssocPion0RecFor);
        GetTriggAssocLists(arrPion0SideMid, listTriggPion0SideMid, listAssocPion0SideMid);
        GetTriggAssocLists(arrPion0SideFor, listTriggPion0SideFor, listAssocPion0SideFor);

        DoCorrelations(arrPion0Mid, listTriggMid, listAssocMid, hCorrMid);
        DoCorrelations(arrPion0For, listTriggFor, listAssocFor, hCorrFor);
        DoCorrelations(arrChargedFor, listTriggChargedFor, listAssocChargedFor, hCorrChargedFor);
        DoCorrelations(arrChargedMid, listTriggChargedMid, listAssocChargedMid, hCorrChargedMid);
        DoCorrelations(arrPion0RecMid, listTriggPion0RecMid, listAssocPion0RecMid, hCorrPion0RecMid);
        DoCorrelations(arrPion0RecFor, listTriggPion0RecFor, listAssocPion0RecFor, hCorrPion0RecFor);
        DoCorrelations(arrPion0SideMid, listTriggPion0SideMid, listAssocPion0SideMid, hCorrPion0SideMid);
        DoCorrelations(arrPion0SideFor, listTriggPion0SideFor, listAssocPion0SideFor, hCorrPion0SideFor);
    }

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    pythia.stat();

    return 0;
}

//-------------------------
//        Functions       |
//-------------------------
void DoCorrelations(TClonesArray *arrPions, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins])
{
    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    for (int it=0; it<nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrPions->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();
        int iTriggBin = GetBin(triggPt, nTriggBins, ptTrigg);
        
        int nAssoc = listAssoc.size();       
        if (nAssoc < 1) continue;

        for (int ia=0; ia<nAssoc; ia++) {
            int iAssoc = listAssoc[ia];

            if (iTrigg==iAssoc) continue; // autocorrelations

            TLorentzVector *lvAssoc = (TLorentzVector*)arrPions->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, nAssocBins, ptAssoc);

            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;

            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta);
        }
    }
}

int GetBin(double arr[], int nArr, double val)
{
    for (int i=0; i<nArr-1; i++) {
        if (arr[i]<=val && val<arr[i+1]) return i;
    }
    return -1;
}

void GetTriggAssocLists(TClonesArray *arrPions, std::vector<int>& listTrigg, std::vector<int>& listAssoc)
{
    int nPions = arrPions->GetEntriesFast();
    if (!nPions) return;

    for (int i=0; i<nPions; i++) {
        TLorentzVector *lvPion = (TLorentzVector*)arrPions->At(i);
        double pt = lvPion->Pt();
        int iTrigg = GetBin(triggPt, nTriggBins, pt);
        int iAssoc = GetBin(assocPt, nAssocBins, pt);
        if (iTrigg >= 0) listTrigg.push_back(i);
        if (iAssoc >= 0) listAssoc.push_back(i);
    }
}

double GetDeltaPhi(double phiTrigg, double phiAssoc)
{
    double dphi = phiTrigg - phiAssoc;
    if (dphi>3.0*TMath::Pi()/2.0) return dphi - 2.0*TMath::Pi();
    if (dphi<-TMath::Pi()/2.0) return dphi + 2.0*TMath::Pi();
    return dphi;
}

// For FoCal:
//      sigma/E = a/sqrt(E) + b = 0.25/sqrt(E) + 0.01
double PhotonEnergySmearing(TRandom3 *rand, double px, double py, double pz)
{
    double a = 0.25;
    double b = 0.01;
    double e = sqrt(px*px + py*py + pz*pz);
    double sigma = sqrt(a*a*e + b*b*e*e);
    double eSmear = rand->Gaus(e, sigma);
    return eSmear;
}

void ReconstructPions(TClonesArray *arrayPhoton, TClonesArray *arrayPi0Rec, TClonesArray *arrayPi0Side)
{
    int nRec = 0, nSide = 0;
    int nPhoton = arrayPhoton->GetEntriesFast();
    for (int i = 0; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector *lv1 = (TLorentzVector*)arrayPhoton->At(i);
            TLorentzVector *lv2 = (TLorentzVector*)arrayPhoton->At(j);
            TLorentzVector lvSum;
            lvSum = ((*lv1) + (*lv2));
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            int ibin = GetBin(pi0Pt, nPi0PtBins, pT);
            if (ibin < 0) continue;
            if (mass > 110. && mass < 170.) new ((*arrayPi0Rec)[nRec++]) TLorentzVector(lvSum);
            if ((mass > 55. && mass < 75.) || (mass > 210. && mass < 280.)) new ((*arrayPi0Side)[nSide++]) TLorentzVector(lvSum);
        }
    }
}

void FillPionMasses(TClonesArray *arrayPhoton, TH1D *hMasses[nPi0PtBins])
{
    int nPhoton = arrayPhoton->GetEntriesFast();
    for (int i = 0; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector *lv1 = (TLorentzVector*)arrayPhoton->At(i);
            TLorentzVector *lv2 = (TLorentzVector*)arrayPhoton->At(j);
            TLorentzVector lvSum;
            lvSum = ((*lv1) + (*lv2));
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            int ibin = GetBin(pi0Pt, nPi0PtBins, pT);
            if (ibin < 0) continue;
            hMasses[ibin]->Fill(mass);
        }
    }   
}
