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
#include "TStopwatch.h"

using namespace Pythia8;

//-------------------------
//         pT bins        |
//-------------------------
const int nTriggBins = 1;
double  triggPt[nTriggBins+1] = {5.0, 15.0};
//double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

const int nAssocBins = 1;
double  assocPt[nAssocBins+1] = {2.0, 5.0};
//double  assocPt[nAssocBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};

// Minimun pT when using leading trigger
double pTmin = 3.0; 

//-------------------------
//  Histogram parameters  |
//-------------------------
const double etaBinWidth = 0.025;
const double phiBinWidth = 0.025;

const double etaTrackerRange = 0.9;
const double etaFocalMin = 3.2;
const double etaFocalMax = 5.8;
const double etaFocalRange = 5.8 - 3.2;
const int nEtaBinTracker = int(etaTrackerRange/etaBinWidth) + 1;
const int nEtaBinFocal = int(etaFocalRange/etaBinWidth) + 1;

const double deltaPhiMin = -TMath::Pi()/2.0;
const double deltaPhiMax = 3.0*TMath::Pi()/2.0;
const int nPhiBin = int((deltaPhiMax-deltaPhiMin)/phiBinWidth) + 1;

const int nIncPtBin = 150;
double logBinsX[nIncPtBin+1], limMin = 0.1, limMax = 100;
const double logBW = (log(limMax) - log(limMin))/nIncPtBin;

const int nIncEtaBin = 150;
const double incEtaRange = 20.0;

const int nPool = 5;

//-------------------------
//        Functions       |
//-------------------------
bool IsTrackerAcceptance(double eta, double etaRange=etaTrackerRange);
bool IsFocalAcceptance(double eta, double etaMin=etaFocalMin, double etaMax=etaFocalMax);
void DoCorrelations(TClonesArray *arrPionsTrigg, TClonesArray *arrPionsAssoc, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins]);
void DoCorrelations(TClonesArray *arrPi0, TLorentzVector *lvTrigg, int itrigg, bool bTriggIncluded, TH2D *hCorr);
int GetBin(double arr[], int nArr, double val);
void GetTriggAssocLists(TClonesArray *arrPi0Trigg, TClonesArray *arrPi0Assoc, std::vector<int>& listTrigg, std::vector<int>& listAssoc);
double GetDeltaPhi(double phiTrigg, double phiAssoc);

double PhotonEnergySmearing(TRandom3 *rand, double px, double py, double pz);TLorentzVector GetPhotonSumVector(TClonesArray *arrayPhoton, int iPhoton1, int iPhoton2);
TLorentzVector GetPhotonSumVector(TClonesArray *arrayPhoton, int iPhoton1, int iPhoton2);
void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0trigg, TClonesArray *arrPi0assoc, bool bMass);
void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0, bool bMass);
void FillPionMasses(TClonesArray *arrPhoton, TH1D *hMassesTrigg[nTriggBins], TH1D *hMassesAssoc[nAssocBins]);
void FillPionMasses(TClonesArray *arrPi0, TH1D *hMassesTrigg, TH1D *hMassesAssoc, int itrigg);
int GetLeadingTriggerIndex(TClonesArray *arrPi0);

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

    TRandom3 *rand = new TRandom3();

    // Basic histograms
    TH1D *hCounter = new TH1D("hCounter", "hCounter", 10, 0, 10);
    
    for (int i = 0; i <= nIncPtBin; i++) logBinsX[i] = limMin*exp(i*logBW);
    
    TH1D *hPionPt = new TH1D("hPionPt", "hPionPt", nIncPtBin, logBinsX); hPionPt->Sumw2();
    TH1D *hPionPtFor = new TH1D("hPionPtFor", "hPionPtFor", nIncPtBin, logBinsX); hPionPtFor->Sumw2();
    TH1D *hPionPtMid = new TH1D("hPionPtMid", "hPionPtMid", nIncPtBin, logBinsX); hPionPtMid->Sumw2();
    
    //TH1D *hChargedHadronPt = new TH1D("hChargedHadronPt", "hChargedHadronPt", nIncPtBin, logBinsX); hChargedHadronPt->Sumw2();
    //TH1D *hChargedHadronPtFor = new TH1D("hChargedHadronPtFor", "hChargedHadronPtFor", nIncPtBin, logBinsX); hChargedHadronPtFor->Sumw2();
    //TH1D *hChargedHadronPtMid = new TH1D("hChargedHadronPtMid", "hChargedHadronPtMid", nIncPtBin, logBinsX); hChargedHadronPtMid->Sumw2();

    TH1D *hPhotonPt = new TH1D("hPhotonPt", "hPhotonPt", nIncPtBin, logBinsX); hPhotonPt->Sumw2();
    TH1D *hPhotonPtFor = new TH1D("hPhotonPtFor", "hPhotonPtFor", nIncPtBin, logBinsX); hPhotonPtFor->Sumw2();
    TH1D *hPhotonPtMid = new TH1D("hPhotonPtMid", "hPhotonPtMid", nIncPtBin, logBinsX); hPhotonPtMid->Sumw2();

    TH1D *hPionEta = new TH1D("hPionEta", "hPionEta", nIncEtaBin, -incEtaRange/2., incEtaRange/2.); hPionEta->Sumw2();
    //TH1D *hChargedHadronEta = new TH1D("hChargedHadronEta", "hChargedHadronEta", nIncEtaBin, -incEtaRange/2., incEtaRange/2.); hChargedHadronEta->Sumw2();

    // Correlation histograms
    TH2D *hCorrMid[nTriggBins][nAssocBins];
    TH2D *hCorrFor[nTriggBins][nAssocBins];
    //TH2D *hCorrChargedMid[nTriggBins][nAssocBins];
    //TH2D *hCorrChargedFor[nTriggBins][nAssocBins];
    
    TH2D *hCorrMassMass[nTriggBins][nAssocBins];
    TH2D *hCorrMassSide[nTriggBins][nAssocBins];
    TH2D *hCorrSideMass[nTriggBins][nAssocBins];
    TH2D *hCorrSideSide[nTriggBins][nAssocBins];
    
    /**TH2D *hCorrMassMassMixed[nTriggBins][nAssocBins];
    TH2D *hCorrMassSideMixed[nTriggBins][nAssocBins];
    TH2D *hCorrSideMassMixed[nTriggBins][nAssocBins];
    TH2D *hCorrSideSideMixed[nTriggBins][nAssocBins];**/

    for (int i = 0; i < nTriggBins; i++) {
        for (int j = 0; j < nAssocBins; j++) {
            hCorrMid[i][j] = new TH2D(Form("hCorrMid%d:%d", i, j), Form("hCorrMid%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -etaTrackerRange/2., etaTrackerRange/2.);
            hCorrMid[i][j]->Sumw2();
            hCorrFor[i][j] = new TH2D(Form("hCorrFor%d:%d", i, j), Form("hCorrFor%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrFor[i][j]->Sumw2();
            /**hCorrChargedMid[i][j] = new TH2D(Form("hCorrChargedMid%d:%d", i, j), Form("hCorrChargedMid%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -etaTrackerRange/2., etaTrackerRange/2.);
            hCorrChargedMid[i][j]->Sumw2();
            hCorrChargedFor[i][j] = new TH2D(Form("hCorrChargedFor%d:%d", i, j), Form("hCorrChargedFor%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrChargedFor[i][j]->Sumw2();**/
            
            hCorrMassMass[i][j] = new TH2D(Form("hCorrMassMass%d:%d", i, j), Form("hCorrMassMass%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrMassMass[i][j]->Sumw2();
            hCorrMassSide[i][j] = new TH2D(Form("hCorrMassSide%d:%d", i, j), Form("hCorrMassSide%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrMassSide[i][j]->Sumw2();
            hCorrSideMass[i][j] = new TH2D(Form("hCorrSideMass%d:%d", i, j), Form("hCorrSideMass%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrSideMass[i][j]->Sumw2();
            hCorrSideSide[i][j] = new TH2D(Form("hCorrSideSide%d:%d", i, j), Form("hCorrSideSide%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrSideSide[i][j]->Sumw2();
           
            /**hCorrMassMassMixed[i][j] = new TH2D(Form("hCorrMassMassMixed%d:%d", i, j), Form("hCorrMassMassMixed%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrMassMassMixed[i][j]->Sumw2();
            hCorrMassSideMixed[i][j] = new TH2D(Form("hCorrMassSideMixed%d:%d", i, j), Form("hCorrMassSideMixed%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrMassSideMixed[i][j]->Sumw2();
            hCorrSideMassMixed[i][j] = new TH2D(Form("hCorrSideMassMixed%d:%d", i, j), Form("hCorrSideMassMixed%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrSideMassMixed[i][j]->Sumw2();
            hCorrSideSideMixed[i][j] = new TH2D(Form("hCorrSideSideMixed%d:%d", i, j), Form("hCorrSideSideMixed%d:%d", i, j), nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange/2., etaFocalRange/2.);
            hCorrSideSideMixed[i][j]->Sumw2();**/
        }
    }

    TH1D *hPi0MassTrigg[nTriggBins];
    TH1D *hPi0SideTrigg[nTriggBins];
    for (int i = 0; i < nTriggBins; i++) {
        hPi0MassTrigg[i] = new TH1D(Form("hPi0MassTrigg%d", i), Form("hPi0MassTrigg%d", i), 301, 0.0, 300.0);
        hPi0SideTrigg[i] = new TH1D(Form("hPi0SideTrigg%d", i), Form("hPi0SideTrigg%d", i), 301, 0.0, 300.0);
    }

    TH1D *hPi0MassAssoc[nAssocBins];
    TH1D *hPi0SideAssoc[nAssocBins];
    for (int i = 0; i < nAssocBins; i++) {
        hPi0MassAssoc[i] = new TH1D(Form("hPi0MassAssoc%d", i), Form("hPi0MassAssoc%d", i), 301, 0.0, 300.0);
        hPi0SideAssoc[i] = new TH1D(Form("hPi0SideAssoc%d", i), Form("hPi0SideAssoc%d", i), 301, 0.0, 300.0);
    }

    // Particle lists
    TClonesArray *arrPion0Mid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0For = new TClonesArray("TLorentzVector", 1500);
    //TClonesArray *arrChargedMid = new TClonesArray("TLorentzVector", 1500);
    //TClonesArray *arrChargedFor = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPhotonMid = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPhotonFor = new TClonesArray("TLorentzVector", 1500);

    TClonesArray *arrPion0MassTrigg = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0SideTrigg = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0MassAssoc = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPion0SideAssoc = new TClonesArray("TLorentzVector", 1500);
    
    // Vectors to store events for event mixing
    //std::vector<vector<int>> assocMassPool, assocSidePool;
    //std::vector<TClonesArray> arrMassPool, arrSidePool;

    //
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        hCounter->Fill(0.5); // Number of events

        //std::cout << "event " << iEvent << std::endl;

        if ( !pythia.next() ) continue;

        int nPion0Mid = 0;
        int nPion0For = 0;
        //int nChargedMid = 0;
        //int nChargedFor = 0;
        int nPhotonMid = 0;
        int nPhotonFor = 0;

        arrPion0Mid->Clear();
        arrPion0For->Clear();
        //arrChargedMid->Clear();
        //arrChargedFor->Clear();
        arrPhotonMid->Clear();
        arrPhotonFor->Clear();
    
        arrPion0MassTrigg->Clear(); 
        arrPion0SideTrigg->Clear();
        arrPion0MassAssoc->Clear(); 
        arrPion0SideAssoc->Clear();

        // Collect particles of interest
        for (int iPart = 0; iPart < pythia.event.size(); iPart++) {

            TLorentzVector lv(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), pythia.event[iPart].e());

            if (pythia.event[iPart].id()==22) { // photons
                double px = pythia.event[iPart].px();
                double py = pythia.event[iPart].py();
                double pz = pythia.event[iPart].pz();
                double e = TMath::Sqrt(px*px + py*py + pz*pz);
                double eSmear = PhotonEnergySmearing(rand, px, py, pz);
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

            if (pythia.event[iPart].id()==111) { // Pions
                double pt = pythia.event[iPart].pT();
                double eta = pythia.event[iPart].eta();

                hPionPt->Fill(pt);
                hPionEta->Fill(eta);
                if (IsTrackerAcceptance(eta)) {
                    hPionPtMid->Fill(pt);
                    new((*arrPion0Mid)[nPion0Mid++]) TLorentzVector(lv);
                }
                if (IsFocalAcceptance(eta)) {
                    hPionPtFor->Fill(pt);
                    new((*arrPion0For)[nPion0For++]) TLorentzVector(lv);
                }
            }

            /**if (pythia.event[iPart].isFinal() && pythia.event[iPart].isHadron() && pythia.event[iPart].isCharged()) { // charged hadrons
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
            }**/
        }
        
        if (nPion0Mid>0) hCounter->Fill(1.5); // number of events with pion0 in mid rapidity
        if (nPion0For>0) hCounter->Fill(2.5); // number of events with pion0 in forward rapidity 

        if (bUseLeading) {
            // Real pions
            int itrigg = GetLeadingTriggerIndex(arrPion0For);
            if (itrigg >= 0) {
                TLorentzVector *lvTriggReal = (TLorentzVector*)arrPion0For->At(itrigg);
                DoCorrelations(arrPion0For, lvTriggReal, itrigg, 1, hCorrFor[0][0]);
            }
            // Reconstructed pions
            ReconstructPions(arrPhotonFor, arrPion0MassAssoc, 1);
            ReconstructPions(arrPhotonFor, arrPion0SideAssoc, 0);
            
            int itriggMass = GetLeadingTriggerIndex(arrPion0MassAssoc);
            int itriggSide = GetLeadingTriggerIndex(arrPion0SideAssoc);  
             
            if (itriggMass > itriggSide) {
                if (itriggMass < 0) continue;
                TLorentzVector *lvTrigg = (TLorentzVector*)arrPion0MassAssoc->At(itriggMass);
                DoCorrelations(arrPion0MassAssoc, lvTrigg, itriggMass, 1, hCorrMassMass[0][0]);
                DoCorrelations(arrPion0SideAssoc, lvTrigg, itriggMass, 0, hCorrMassSide[0][0]);
            } else {
                if (itriggSide < 0) continue;
                TLorentzVector *lvTrigg = (TLorentzVector*)arrPion0SideAssoc->At(itriggSide);
                DoCorrelations(arrPion0MassAssoc, lvTrigg, itriggSide, 0, hCorrSideMass[0][0]);
                DoCorrelations(arrPion0SideAssoc, lvTrigg, itriggSide, 1, hCorrSideSide[0][0]);
            }  
        } else {
            // Real pions
            std::vector<int> listTriggMid, listTriggFor, listTriggChargedMid, listTriggChargedFor;
            std::vector<int> listAssocMid, listAssocFor, listAssocChargedMid, listAssocChargedFor;
        
            GetTriggAssocLists(arrPion0Mid, arrPion0Mid, listTriggMid, listAssocMid);
            GetTriggAssocLists(arrPion0For, arrPion0For, listTriggFor, listAssocFor);
        
            if (listTriggFor.size()>0) hCounter->Fill(3.5); // number of events with trigger pion0 in forward rapidity 
            if (listTriggFor.size()>0) hCounter->Fill(4.5, (int)listTriggFor.size()); // number of real trigger pion0 in forward rapidity 
            if (listTriggFor.size() && listAssocFor.size()>0) hCounter->Fill(5.5, (int)listAssocFor.size()); // number of real assoc pion0 in forward rapidity 
        
            DoCorrelations(arrPion0Mid, arrPion0Mid, listTriggMid, listAssocMid, hCorrMid);
            DoCorrelations(arrPion0For, arrPion0For, listTriggFor, listAssocFor, hCorrFor);
        
            // Reconstructed pions       
            std::vector<int> listTriggMass, listTriggSide, listAssocMass, listAssocSide;

            ReconstructPions(arrPhotonFor, arrPion0MassTrigg, arrPion0MassAssoc, 1);
            ReconstructPions(arrPhotonFor, arrPion0SideTrigg, arrPion0SideAssoc, 0);

            int nTriggMass = arrPion0MassTrigg->GetEntriesFast();
            int nTriggSide = arrPion0SideTrigg->GetEntriesFast();

            if (nTriggMass>0) FillPionMasses(arrPhotonFor, hPi0MassTrigg, hPi0MassAssoc);
            if (nTriggSide>0) FillPionMasses(arrPhotonFor, hPi0SideTrigg, hPi0SideAssoc);

            GetTriggAssocLists(arrPion0MassTrigg, arrPion0MassAssoc, listTriggMass, listAssocMass);
            GetTriggAssocLists(arrPion0SideTrigg, arrPion0SideAssoc, listTriggSide, listAssocSide);
               
            DoCorrelations(arrPion0MassTrigg, arrPion0MassAssoc, listTriggMass, listAssocMass, hCorrMassMass);
            DoCorrelations(arrPion0MassTrigg, arrPion0SideAssoc, listTriggMass, listAssocSide, hCorrMassSide);
            DoCorrelations(arrPion0SideTrigg, arrPion0MassAssoc, listTriggSide, listAssocMass, hCorrSideMass);
            DoCorrelations(arrPion0SideTrigg, arrPion0SideAssoc, listTriggSide, listAssocSide, hCorrSideSide);
        }
        // Event mixing
        /**if (iEvent >= nPool) {
           for (int ipool = 0; ipool < nPool; ipool++) {
               DoCorrelations(arrPion0MassTrigg, &arrMassPool[ipool], listTriggMass, assocMassPool[ipool], hCorrMassMassMixed);
               DoCorrelations(arrPion0MassTrigg, &arrSidePool[ipool], listTriggMass, assocSidePool[ipool], hCorrMassSideMixed);
               DoCorrelations(arrPion0SideTrigg, &arrMassPool[ipool], listTriggSide, assocMassPool[ipool], hCorrSideMassMixed);
               DoCorrelations(arrPion0SideTrigg, &arrSidePool[ipool], listTriggSid, assocSidePool[ipool], hCorrSideSideMixed);
           }

           assocMassPool.erase(assocMassPool.begin());
           assocSidePool.erase(assocSidePool.begin());
           arrMassPool.erase(arrMassPool.begin());
           arrSidePool.erase(arrSidePool.begin());
        }
                
        assocMassPool.push_back(listAssocMass);
        assocSidePool.push_back(listAssocSide);
        arrMassPool.push_back(*arrPion0MassAssoc);
        arrSidePool.push_back(*arrPion0SideAssoc);**/

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

void DoCorrelations(TClonesArray *arrPionsTrigg, TClonesArray *arrPionsAssoc, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[nTriggBins][nAssocBins])
{
    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    for (int it=0; it<nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrPionsTrigg->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();
        int iTriggBin = GetBin(triggPt, nTriggBins, ptTrigg);
        
        int nAssoc = listAssoc.size();       
        if (nAssoc < 1) continue;

        for (int ia=0; ia<nAssoc; ia++) {
            int iAssoc = listAssoc[ia];

            if (iTrigg==iAssoc) continue; // autocorrelations

            TLorentzVector *lvAssoc = (TLorentzVector*)arrPionsAssoc->At(iAssoc);
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

void DoCorrelations(TClonesArray *arrPi0, TLorentzVector *lvTrigg, int itrigg, bool bTriggIncluded, TH2D *hCorr)
{
    int nPi0 = arrPi0->GetEntriesFast();
    if (!nPi0) return;

    double ptTrigg = lvTrigg->Pt();
    double phiTrigg = lvTrigg->Phi();
    double etaTrigg = lvTrigg->Eta();
    
    for (int i = 0; i < nPi0; i++) {
        if (bTriggIncluded && i == itrigg) continue;
        TLorentzVector *lvAssoc = (TLorentzVector*)arrPi0->At(i);
        double ptAssoc = lvAssoc->Pt();
        double phiAssoc = lvAssoc->Phi();
        double etaAssoc = lvAssoc->Eta();

        double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
        double deta = etaTrigg - etaAssoc;

        hCorr->Fill(dphi, deta);
    }
}

int GetBin(double arr[], int nArr, double val)
{
    for (int i=0; i<nArr; i++) {
        if (arr[i]<=val && val<arr[i+1]) return i;
    }
    return -1;
}

void GetTriggAssocLists(TClonesArray *arrPi0Trigg, TClonesArray *arrPi0Assoc, std::vector<int>& listTrigg, std::vector<int>& listAssoc)
{
    int nPi0Trigg = arrPi0Trigg->GetEntriesFast();
    int nPi0Assoc = arrPi0Assoc->GetEntriesFast();
    
    //if (!nPi0Trigg) return;
    //if (!nPi0Assoc) return;

    for (int i=0; i<nPi0Trigg; i++) {
        TLorentzVector *lvPion = (TLorentzVector*)arrPi0Trigg->At(i);
        double pt = lvPion->Pt();
        int iTrigg = GetBin(triggPt, nTriggBins, pt);
        if (iTrigg >= 0) listTrigg.push_back(i);
    }
   
    for (int i=0; i<nPi0Assoc; i++) {
        TLorentzVector *lvPion = (TLorentzVector*)arrPi0Assoc->At(i);
        double pt = lvPion->Pt();
        int iAssoc = GetBin(assocPt, nAssocBins, pt);
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
//      sigma/E = a/sqrt(E) + b = 0.27/sqrt(E) + 0.01
double PhotonEnergySmearing(TRandom3 *rand, double px, double py, double pz)
{
    const double eCut = 0.01;
    double e = sqrt(px*px + py*py + pz*pz);
    if (e < eCut) return e;

    double a = 0.27;
    double b = 0.01;
    double sigma = sqrt(a*a*e + b*b*e*e);
    double eSmear = rand->Gaus(e, sigma);
    return eSmear;
}

TLorentzVector GetPhotonSumVector(TClonesArray *arrPhoton, int iPhoton1, int iPhoton2)
{
    TLorentzVector *lv1 = (TLorentzVector*)arrPhoton->At(iPhoton1);
    TLorentzVector *lv2 = (TLorentzVector*)arrPhoton->At(iPhoton2);
    TLorentzVector lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}

void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0trigg, TClonesArray *arrPi0assoc, bool bMass)
{
    int nTrigg = 0, nAssoc = 0;
    int nPhoton = arrPhoton->GetEntriesFast();

    for (int i = 1; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector lvSum = GetPhotonSumVector(arrPhoton, i, j);
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            int iTriggBin = GetBin(triggPt, nTriggBins, pT);
            int iAssocBin = GetBin(assocPt, nAssocBins, pT);
            bool bIsInWindow = bMass ? IsMassWindow(mass) : IsSideband(mass); 
            if (bIsInWindow) {
                if (iTriggBin >= 0) new ((*arrPi0trigg)[nTrigg++]) TLorentzVector(lvSum);
                if (iAssocBin >= 0) new ((*arrPi0assoc)[nAssoc++]) TLorentzVector(lvSum);
            }
        }
    }
}

void ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0, bool bMass)
{
    int nPi0 = 0;
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector lvSum = GetPhotonSumVector(arrPhoton, i, j);
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            if (pT<pTmin) continue;
            bool bIsInWindow = bMass ? IsMassWindow(mass) : IsSideband(mass); 
            if (bIsInWindow) new ((*arrPi0)[nPi0++]) TLorentzVector(lvSum);
        }
    }
}

void FillPionMasses(TClonesArray *arrPhoton, TH1D *hMassesTrigg[nTriggBins], TH1D *hMassesAssoc[nAssocBins])
{
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector lvSum = GetPhotonSumVector(arrPhoton, i, j);
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            int iTriggBin = GetBin(triggPt, nTriggBins, pT);
            int iAssocBin = GetBin(assocPt, nAssocBins, pT);
            if (iTriggBin >= 0) hMassesTrigg[iTriggBin]->Fill(mass);
            if (iAssocBin >= 0) hMassesAssoc[iAssocBin]->Fill(mass);
        }
    }   
}

void FillPionMasses(TClonesArray *arrPi0, TH1D *hMassesTrigg, TH1D *hMassesAssoc, int itrigg)
{
    int nPi0 = arrPi0->GetEntriesFast();
    for (int i = 0; i < nPi0; i++) {
        TLorentzVector *lvPi0 = (TLorentzVector*)arrPi0->At(i);
        double mass = 1000.*lvPi0->M();
        if (i == itrigg) {
            hMassesTrigg->Fill(mass);
        } else {
            hMassesAssoc->Fill(mass);
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
        if (pTmax < pT) {
            pTmax = pT;
            itrigg = i;
        }
    }
    return itrigg;
}
