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

#include "include/AliJBaseTrack.h"
#include "include/AliJHMRConst.h"
#include "include/AliJHMRHist.h"
#include "include/AliJHMRCorr.h"
#include "include/AliJHMRPythiaCatalyst.h"

using namespace Pythia8;

int main(int argc, char *argv[]) {

    if (argc==1) {
        cout << "Usage : ./focalPionCorrelation output.root bUseLeading bDebugOn pythiaSettings.cmnd poolsize seed" << endl;
        return 0;
    }

    TStopwatch timer;
    timer.Start();

    TString outFileName = argc > 1 ? argv[1] : "output.root";
    int bUseLeading = argc > 2 ? atol(argv[2]) : 0;
    int bDebugOn = argc > 3 ? atol(argv[3]) : 0;
    TString pythiaSettings = argc > 4 ? argv[4] : "PythiaHard.cmnd";
    int poolsize = argc > 5 ? atol(argv[5]) : 10;
    int seed = argc > 6 ? atol(argv[6]) : 0;

    detector det = kJFoCal;

    TFile *fOut = new TFile(outFileName, "RECREATE");

    Pythia pythia;

    // Initialise pythia
    pythia.readFile(pythiaSettings.Data());
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.init();

    int nEvents = pythia.mode("Main:numberOfEvents");

    AliJHMRHist *fHistos = new AliJHMRHist();
    fHistos->CreateHistos(fOut, det);

    AliJHMRPythiaCatalyst *fCatalyst = new AliJHMRPythiaCatalyst(pythia.event, fHistos);
    AliJHMRCorr *fCorr = new AliJHMRCorr(fHistos, det);

    TClonesArray *arrPhotonFor = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Real = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Peak = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Side = new TClonesArray("AliJBaseTrack", 1500);

    std::vector<TClonesArray*> arrPi0PeakMixed(poolsize);
    std::vector<TClonesArray*> arrPi0SideMixed(poolsize);

    for (int ipool = 0; ipool < poolsize; ipool++) {
        arrPi0PeakMixed[ipool] = new TClonesArray("AliJBaseTrack", 1500);
        arrPi0SideMixed[ipool] = new TClonesArray("AliJBaseTrack", 1500);
    }

    fOut->cd();

    //
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        fHistos->hCounter->Fill(0.5); // number of events

        if (bDebugOn)
            std::cout << "\nEvent " << iEvent << std::endl;

        if ( !pythia.next() ) continue;

        fCatalyst->InitializeEvent();
        fCatalyst->GetParticles(det);
        arrPhotonFor = fCatalyst->GetParticleList(kJDecayPhoton);
        arrPi0Real = fCatalyst->GetParticleList(kJPi0);

        fCorr->SmearEnergies(arrPhotonFor);

        fHistos->FillPtEta(kJDecayPhoton, arrPhotonFor);
        fHistos->FillPtEta(kJPi0, arrPi0Real);

        fCorr->ReconstructPions(arrPhotonFor, arrPi0Peak, 1);
        fCorr->ReconstructPions(arrPhotonFor, arrPi0Side, 0);

        if (bDebugOn)
            std::cout << "Number of Pi0 (real) : " << arrPi0Real->GetEntriesFast()
                      << "\t (rec, peak) : " << arrPi0Peak->GetEntriesFast()
                      << "\t gamma : " << arrPhotonFor->GetEntriesFast() << std::endl;

        std::vector<int> listTriggReal, listAssocReal, listTriggPeak, listTriggSide, listAssocPeak, listAssocSide;
        int binsWithTriggReal[NTRIGGBINS+1] = {0}, binsWithTriggPeak[NTRIGGBINS+1] = {0}, binsWithTriggSide[NTRIGGBINS+1] = {0};
        fCorr->GetTriggAssocLists(arrPi0Real, listTriggReal, listAssocReal, binsWithTriggReal, bUseLeading);
        fCorr->GetTriggAssocLists(arrPi0Peak, listTriggPeak, listAssocPeak, binsWithTriggPeak, bUseLeading);
        fCorr->GetTriggAssocLists(arrPi0Side, listTriggSide, listAssocSide, binsWithTriggSide, bUseLeading);

        fCorr->FillPionMasses(arrPhotonFor, binsWithTriggPeak, binsWithTriggSide);
        fCorr->FillRealTriggers(arrPi0Real, listTriggReal);

        fCorr->DoCorrelations(arrPi0Real, listTriggReal, listAssocReal, fHistos->hCorrFor, 0);
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

        // Construct & save true correlation components f_SS, f_SB, f_BS, f_BB
        fCorr->ConstructTrueCorrComponents(arrPi0Peak, listTriggPeak, listAssocPeak, 0);

        // Mixed event : take triggers from this event, associated from previous
        if (iEvent >= poolsize) {
            for (int ipool = 0; ipool < poolsize; ipool++) {
                std::vector<int> listTriggPeakMixed, listTriggSideMixed, listAssocPeakMixed, listAssocSideMixed;
                int binsWithTriggPeakMixed[NTRIGGBINS+1] = {0}, binsWithTriggSideMixed[NTRIGGBINS+1] = {0};
                fCorr->GetTriggAssocLists(arrPi0PeakMixed[ipool], listTriggPeakMixed, listAssocPeakMixed, binsWithTriggPeakMixed, bUseLeading);
                fCorr->GetTriggAssocLists(arrPi0SideMixed[ipool], listTriggSideMixed, listAssocSideMixed, binsWithTriggSideMixed, bUseLeading);
                fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0PeakMixed[ipool], listAssocPeakMixed, fHistos->hCorrMassMassMixed, 0, 0);
                fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0SideMixed[ipool], listAssocSideMixed, fHistos->hCorrSideSideMixed, 0, 0);
                fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0SideMixed[ipool], listAssocSideMixed, fHistos->hCorrMassSideMixed, 0, 0);
                fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0PeakMixed[ipool], listAssocPeakMixed, fHistos->hCorrSideMassMixed, 0, 0);

            }
            // Remove the first from the pool and add new array to the pool
            arrPi0PeakMixed.erase(arrPi0PeakMixed.begin());
            arrPi0PeakMixed.push_back((TClonesArray*)arrPi0Peak->Clone());
            arrPi0SideMixed.erase(arrPi0SideMixed.begin());
            arrPi0SideMixed.push_back((TClonesArray*)arrPi0Side->Clone());
        } else { // Create pools"
            arrPi0PeakMixed[iEvent] = (TClonesArray*)arrPi0Peak->Clone();
            arrPi0SideMixed[iEvent] = (TClonesArray*)arrPi0Side->Clone();
        }

        arrPhotonFor->Clear("C");
        arrPi0Real->Clear("C");
        arrPi0Peak->Clear("C");
        arrPi0Side->Clear("C");
    }

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    pythia.stat();

    timer.Print();

    return 0;
}
