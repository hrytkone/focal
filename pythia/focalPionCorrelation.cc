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
#include "TObjectTable.h"

#include "include/AliJBaseTrack.h"
#include "include/AliJHMRConst.h"
#include "include/AliJHMRHist.h"
#include "include/AliJHMRCorr.h"
#include "include/AliJHMRPythiaCatalyst.h"
#include "include/AliJHMRGeantCatalyst.h"

using namespace Pythia8;

int main(int argc, char *argv[]) {

    if (argc==1) {
        cout << "Usage : ./focalPionCorrelation output.root bUseLeading bDebugOn pythiaSettings.cmnd poolsize seed bUseSim siminput" << endl;
        return 0;
    }

    TStopwatch timer;
    timer.Start();

    TString outFileName    = argc > 1 ? argv[1] : "output.root";
    int bUseLeading        = argc > 2 ? atol(argv[2]) : 0;
    int bDebugOn           = argc > 3 ? atol(argv[3]) : 0;
    TString pythiaSettings = argc > 4 ? argv[4] : "PythiaHard.cmnd";
    int poolsize           = argc > 5 ? atol(argv[5]) : 10;
    int seed               = argc > 6 ? atol(argv[6]) : 0;
    Int_t bUseSim          = argc > 7 ? atol(argv[7]) : 0;
    TString siminput       = argc > 8 ? argv[8] : "input.root";

    std::cout << "Simulation parameters : " << std::endl;
    std::cout << "\tOutput         : \t" << outFileName << std::endl;
    std::cout << "\tUse leading    : \t" << bUseLeading << std::endl;
    std::cout << "\tDebug          : \t" << bDebugOn << std::endl;
    std::cout << "\tPythia config  : \t" << pythiaSettings << std::endl;
    std::cout << "\tPool size      : \t" << poolsize << std::endl;
    std::cout << "\tSeed           : \t" << seed << std::endl;
    std::cout << "\tUse sim input  : \t" << bUseSim << std::endl;
    std::cout << "\tSim input file : \t" << siminput << std::endl;
    
    detector det = kJFoCal;

    TFile *fOut = new TFile(outFileName, "RECREATE");

    Pythia pythia;

    AliJHMRHist *fHistos = new AliJHMRHist();
    fHistos->CreateHistos(fOut, det);

    AliJHMRPythiaCatalyst *fCatalyst = new AliJHMRPythiaCatalyst(pythia.event, fHistos);
    AliJHMRGeantCatalyst *fCatalystG = new AliJHMRGeantCatalyst(siminput, fHistos);

    AliJHMRCorr *fCorr = new AliJHMRCorr(fHistos, det);

    TClonesArray *arrPhotonFor = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Real = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Peak = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Side = new TClonesArray("AliJBaseTrack", 1500);

    TClonesArray* arrPi0PeakMixed[poolsize];
    TClonesArray* arrPi0SideMixed[poolsize];

    for (int ipool = 0; ipool < poolsize; ipool++) {
        arrPi0PeakMixed[ipool] = new TClonesArray("AliJBaseTrack", 1500);
        arrPi0SideMixed[ipool] = new TClonesArray("AliJBaseTrack", 1500);
    }

    int nEvents = 0;
    if (bUseSim) {
        if (!fCatalystG->LoadInput()) return 1;
        nEvents = fCatalystG->GetNumberOfEvents();
    } else {
        // Initialise pythia
        pythia.readFile(pythiaSettings.Data());
        pythia.readString("Random:setSeed = on");
        pythia.readString(Form("Random:seed = %d", seed));
        pythia.init();
        nEvents = pythia.mode("Main:numberOfEvents");
    }

    std::cout << "Number of events : " << nEvents << std::endl;

    //
    // Loop over events
    //
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        fHistos->hCounter->Fill(0.5); // number of events

        if (bDebugOn)
            std::cout << "\nEvent " << iEvent << std::endl;

        if (bUseSim) { // Get stuff from detector simulation file
            fCatalystG->GetEvent(iEvent);
            arrPhotonFor = fCatalystG->GetClusters();
            arrPi0Real   = fCatalystG->GetPi0True();
        } else { // Get stuff from pythia
            if ( !pythia.next() ) continue;
            fCatalyst->InitializeEvent();
            fCatalyst->GetParticles(det);
            arrPhotonFor = fCatalyst->GetParticleList(kJDecayPhoton);
            arrPi0Real   = fCatalyst->GetParticleList(kJPi0);
            fCorr->SmearEnergies(arrPhotonFor);
        }

        int nTrueFromPeak = fCorr->ReconstructPions(arrPhotonFor, arrPi0Peak, det, 1);
        fCorr->ReconstructPions(arrPhotonFor, arrPi0Side, det, 0);

        fHistos->FillPtEta(kJDecayPhoton, arrPhotonFor);
        fHistos->FillPtEta(kJPi0, arrPi0Real);
        fHistos->FillPtEta(kJRecPi0, arrPi0Peak);
        
        if (bDebugOn)
            std::cout << "Number of Pi0 (real) : " << arrPi0Real->GetEntriesFast()
                      << "\t (rec, peak) : " << arrPi0Peak->GetEntriesFast()
                      << "\t (rec, peak) true : " << nTrueFromPeak
                      << "\t gamma : " << arrPhotonFor->GetEntriesFast() << std::endl;

        std::vector<int> listTriggReal, listAssocReal, listTriggPeak, listTriggSide, listAssocPeak, listAssocSide;
        int binsWithTriggReal[NTRIGGBINS+1] = {0}, binsWithTriggPeak[NTRIGGBINS+1] = {0}, binsWithTriggSide[NTRIGGBINS+1] = {0};
        fCorr->GetTriggAssocLists(arrPi0Real, listTriggReal, listAssocReal, binsWithTriggReal, bUseLeading);
        fCorr->GetTriggAssocLists(arrPi0Peak, listTriggPeak, listAssocPeak, binsWithTriggPeak, bUseLeading);
        fCorr->GetTriggAssocLists(arrPi0Side, listTriggSide, listAssocSide, binsWithTriggSide, bUseLeading);

        fCorr->FillPionMasses(arrPhotonFor, binsWithTriggPeak, binsWithTriggSide, det);
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
                fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0PeakMixed[ipool], listAssocPeakMixed, fHistos->hCorrMassMassMixed, 1, 1);
                fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0SideMixed[ipool], listAssocSideMixed, fHistos->hCorrSideSideMixed, 0, 0);
                fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0SideMixed[ipool], listAssocSideMixed, fHistos->hCorrMassSideMixed, 1, 0);
                fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0PeakMixed[ipool], listAssocPeakMixed, fHistos->hCorrSideMassMixed, 0, 1);

            }
            // Remove the first from the pool and add new array to the pool
            for (int ipool = 0; ipool < poolsize-1; ipool++) {
                *arrPi0PeakMixed[ipool] = *arrPi0PeakMixed[ipool+1];
                *arrPi0SideMixed[ipool] = *arrPi0SideMixed[ipool+1];
            }
            *arrPi0PeakMixed[poolsize-1] = *arrPi0Peak;
            *arrPi0SideMixed[poolsize-1] = *arrPi0Side;
        } else { // Create pools
            *arrPi0PeakMixed[iEvent] = *arrPi0Peak;
            *arrPi0SideMixed[iEvent] = *arrPi0Side;
        }

        arrPhotonFor->Clear("C");
        arrPi0Real->Clear("C");
        arrPi0Peak->Clear("C");
        arrPi0Side->Clear("C");
    }

    fOut->cd();
    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    if (!bUseSim) pythia.stat();

    timer.Print();

    return 0;
}
