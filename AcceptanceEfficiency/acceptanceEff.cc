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

#define NINCPTBIN 150
#define NDET 2

using namespace Pythia8;

const double detEta[NDET][2] = {
    {3.2, 5.8},
	{2.6, 4.0}
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

    TH1D *hPionPtFor = new TH1D("hPionPtFor", "hPionPtFor", NINCPTBIN, logBinsX); hPionPtFor->Sumw2();
    TH1D *hPionPtDetected = new TH1D("hPionPtDetected", "hPionPtDetected", NINCPTBIN, logBinsX); hPionPtDetected->Sumw2();
    TH1D *hPionPtRatio;// = new TH1D("hPionPtRatio", "hPionPtRatio", NINCPTBIN, logBinsX); hPionPtRatio->Sumw2();

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
                hPionPtFor->Fill(pythia.event[partIdx].pT());

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
                        if ( pythia.event[idDaughter1].eta() > detEta[idet][0] && pythia.event[idDaughter2].eta() < detEta[idet][1] )
                            hPionPtDetected->Fill(pythia.event[partIdx].pT());
                    }
                    //else {
                    //    std::cout << "Daughters no photons but " << pythia.event[idDaughter1].id() << " and " << pythia.event[idDaughter2].id() << ", skip" << std::endl;
                    //}
                }
            }
        }
    }

    hPionPtRatio = (TH1D*)hPionPtDetected->Clone("hPionPtRatio");
    hPionPtRatio->Divide(hPionPtFor);

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    pythia.stat();

    timer.Print();

    return 0;
}
