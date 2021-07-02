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

#include "include/AliJHMRConst.h"
#include "include/AliJHMRHist.h"
#include "include/AliJHMRCorr.h"
#include "include/AliJHMRPythiaCatalyst.h"

using namespace Pythia8;

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
    fHistos->CreateHistos(fOut);

    AliJHMRPythiaCatalyst *fCatalyst = new AliJHMRPythiaCatalyst(pythia.event, fHistos);
    AliJHMRCorr *fCorr = new AliJHMRCorr();

    TClonesArray *arrPhotonFor = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPi0Peak = new TClonesArray("TLorentzVector", 1500);
    TClonesArray *arrPi0Side = new TClonesArray("TLorentzVector", 1500);

    fOut->cd();

    // 
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        //std::cout << "event " << iEvent << std::endl;

        if ( !pythia.next() ) continue;
        
        fCatalyst->InitializeEvent();
        fCatalyst->GetParticles(kJDecayPhoton);
        arrPhotonFor = fCatalyst->GetInputList();

        fCorr->ReconstructPions(arrPhotonFor, arrPi0Peak, 1);
        fCorr->ReconstructPions(arrPhotonFor, arrPi0Side, 0);

        std::vector<int> listTriggReal, listAssocReal, listTriggPeak, listTriggSide, listAssocPeak, listAssocSide;
        int binsWithTriggReal[NTRIGGBINS+1] = {0}, binsWithTriggPeak[NTRIGGBINS+1] = {0}, binsWithTriggSide[NTRIGGBINS+1] = {0};
        //fCorr->GetTriggAssocLists(arrPi0For, listTriggReal, listAssocReal, binsWithTriggReal, bUseLeading); 
        fCorr->GetTriggAssocLists(arrPi0Peak, listTriggPeak, listAssocPeak, binsWithTriggPeak, bUseLeading); 
        fCorr->GetTriggAssocLists(arrPi0Side, listTriggSide, listAssocSide, binsWithTriggSide, bUseLeading); 
    
        fCorr->FillPionMasses(arrPhotonFor, fHistos, binsWithTriggPeak, binsWithTriggSide);
        //fCorr->FillRealTriggers(fHistos, arrPi0For, listTriggReal);
       
        //fCorr->DoCorrelations(arrPi0For, listTriggReal, listAssocReal, hCorrFor, 0);
        fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, listAssocPeak, fHistos->hCorrMassMass, 1);
        fCorr->DoCorrelations(arrPi0Side, listTriggSide, listAssocSide, fHistos->hCorrSideSide, 0);
        if (bUseLeading) {
            int isPeakTriggLarger = fCorr->GetLargerTrigg(arrPi0Peak, listTriggPeak, arrPi0Side, listTriggSide);
            if (isPeakTriggLarger) {
                fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, fHistos->hCorrMassSide, 1, 0);
            } else {
                fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, fHistos->hCorrSideMass, 0, 1);
            }
        } else {
            fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, fHistos->hCorrMassSide, 1, 0);
            fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, fHistos->hCorrSideMass, 0, 1);
        }

        arrPhotonFor->Clear("C");
        arrPi0Peak->Clear("C");
        arrPi0Side->Clear("C");
    }

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    pythia.stat();

    timer.Print();

    return 0;
}