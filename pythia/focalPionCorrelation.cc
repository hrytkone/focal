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
#include "AliJHMRConst.h"

using namespace Pythia8;

//-------------------------
//        Functions       |
//-------------------------
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

    AliJHMRHist *fHistos = new AliJHMRHist();
    AliJHMRPythiaCatalyst *fCatalyst = new AliJHMRPythiaCatalyst(histos);
    AliJHMRCorr *fCorr = new AliJHMRCorr(histos);

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
        
        TClonesArray *arrPhotonFor = fCatalyst->GetParticles(kJDecayPhoton);

        fCorr->ReconstructPions(arrPhotonFor, arrPi0Peak, 1);
        fCorr->ReconstructPions(arrPhotonFor, arrPi0Side, 0);

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
