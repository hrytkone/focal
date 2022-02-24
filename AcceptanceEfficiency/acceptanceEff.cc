#include <vector>
#include <cstring>

#include "Pythia8/Pythia.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TStopwatch.h"

//#define NINCPTBIN 150
#define NINCPTBIN 100
#define NETABIN 13
#define NDET 2

int GetEtaBin(double eta, int idet);

using namespace Pythia8;

const double detEta[NDET][2] = {
    {3.2, 5.8}, // FoCal
    {2.6, 4.0}  // STAR
};

int main(int argc, char *argv[]) {

    if (argc==1) {
        cout << "Usage : ./acceptanceEff idet output.root pythiaSettings.cmnd seed" << endl;
        return 0;
    }

    TStopwatch timer;
    timer.Start();

    int idet = argc > 1 ? atol(argv[1]) : 0;
    TString outFileName = argc > 2 ? argv[2] : "output.root";
    TString pythiaSettings = argc > 3 ? argv[3] : "PythiaHard.cmnd";
    int seed = argc > 4 ? atol(argv[4]) : 0;

    TFile *fOut = new TFile(outFileName, "RECREATE");

    double logBinsX[NINCPTBIN+1], limMin = 0.1, limMax = 100;
    double logBW = (log(limMax) - log(limMin))/NINCPTBIN;
    for (int i = 0; i <= NINCPTBIN; i++) logBinsX[i] = limMin*exp(i*logBW);

    double etaRange = detEta[idet][1] - detEta[idet][0];
    double etaBinWidth = 0.025;
    int nEtaBin = int(etaRange/etaBinWidth) + 1;

    TH1D *hPionPtFor[NETABIN];
    TH1D *hPionPtDetected[NETABIN];
    TH1D *hPionPtRatio[NETABIN];
    TH2D *hPionEtaPtFor;
    TH2D *hPionEtaPtDetected;
    TH2D *hPionEtaPtRatio;

    for (int ieta=0; ieta<NETABIN; ieta++) {
        hPionPtFor[ieta] = new TH1D(Form("hPionPtFor_%d", ieta), "hPionPtFor", NINCPTBIN, logBinsX); hPionPtFor[ieta]->Sumw2();
        hPionPtDetected[ieta] = new TH1D(Form("hPionPtDetected_%d", ieta), "hPionPtDetected", NINCPTBIN, logBinsX); hPionPtDetected[ieta]->Sumw2();
    }
    hPionEtaPtFor = new TH2D("hPionEtaPtFor", "hPionEtaPtFor", nEtaBin, detEta[idet][0], detEta[idet][1], NINCPTBIN, logBinsX);
    hPionEtaPtDetected = new TH2D("hPionEtaPtForDetected", "hPionEtaPtForDetected", nEtaBin, detEta[idet][0], detEta[idet][1], NINCPTBIN, logBinsX);

    Pythia pythia;

    // Initialise pythia
    pythia.readFile(pythiaSettings.Data());
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.init();

    int nEvents = pythia.mode("Main:numberOfEvents");

    fOut->cd();

    //
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        if ( !pythia.next() ) continue;

	    for (int partIdx = 0; partIdx < pythia.event.size(); partIdx++) {
            double trEta = pythia.event[partIdx].eta();
    		if ( trEta < detEta[idet][0] || trEta > detEta[idet][1] ) continue;

    		if ( pythia.event[partIdx].id() == 111 ) {
                int ibin = GetEtaBin(pythia.event[partIdx].eta(), idet);
                hPionPtFor[ibin]->Fill(pythia.event[partIdx].pT());
                hPionEtaPtFor->Fill(pythia.event[partIdx].eta(), pythia.event[partIdx].pT());

                int idDaughter1 = pythia.event[partIdx].daughter1();
                int idDaughter2 = pythia.event[partIdx].daughter2();

                if ( idDaughter1==0 && idDaughter2==0 ) {
                    std::cout << "No decay, skip" << std::endl;
                    continue;
                } else if ( idDaughter1==idDaughter2 ) {
                    std::cout << "No decay, just change of momentum, skip" << std::endl;
                    continue;
                } else if ( idDaughter1>0 && idDaughter2==0 ) {
                    std::cout << "Only one daughter, skip" << std::endl;
                    continue;
                } else {
                    if ( pythia.event[idDaughter1].id()==22 && pythia.event[idDaughter2].id()==22 ) {
                        if ( pythia.event[idDaughter1].eta() > detEta[idet][0] && pythia.event[idDaughter2].eta() < detEta[idet][1] ) {
                            hPionPtDetected[ibin]->Fill(pythia.event[partIdx].pT());
                            hPionEtaPtDetected->Fill(pythia.event[partIdx].eta(), pythia.event[partIdx].pT());
                        }
                    }
                    //else {
                    //    std::cout << "Daughters no photons but " << pythia.event[idDaughter1].id() << " and " << pythia.event[idDaughter2].id() << ", skip" << std::endl;
                    //}
                }
            }
        }
    }

    for (int ieta=0; ieta<NETABIN; ieta++) {
        hPionPtRatio[ieta] = (TH1D*)hPionPtDetected[ieta]->Clone(Form("hPionPtRatio_%d", ieta));
        hPionPtRatio[ieta]->Divide(hPionPtFor[ieta]);
    }

    hPionEtaPtRatio = (TH2D*)hPionEtaPtDetected->Clone("hPionEtaPtRatio");
    hPionEtaPtRatio->Divide(hPionEtaPtFor);

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    pythia.stat();

    timer.Print();

    return 0;
}

int GetEtaBin(double eta, int idet)
{
    double step = (detEta[idet][1]-detEta[idet][0])/NETABIN;
    double binEdge = detEta[idet][0];
    for (int i=0; i<NETABIN; i++) {
        if (eta>binEdge && eta<=binEdge+step) return i;
        binEdge += step;
    }
    return -1;
}
