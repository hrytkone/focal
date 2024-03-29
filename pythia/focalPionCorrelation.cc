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
        cout << "Usage : ./focalPionCorrelation <output.root> <bUseLeading> <bDebugOn> <pythiaSettings.cmnd> <poolsize> <seed> <bUseSim> <siminput> <evStart> <evEnd>" << endl;
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
    int evStart            = argc > 9 ? atol(argv[9]) : 0;
    int nEvents            = argc > 10 ? atol(argv[10]) : -1;

    std::cout << "Simulation parameters : " << std::endl;
    std::cout << "\tOutput         : \t" << outFileName << std::endl;
    std::cout << "\tUse leading    : \t" << bUseLeading << std::endl;
    std::cout << "\tDebug          : \t" << bDebugOn << std::endl;
    std::cout << "\tPythia config  : \t" << pythiaSettings << std::endl;
    std::cout << "\tPool size      : \t" << poolsize << std::endl;
    std::cout << "\tSeed           : \t" << seed << std::endl;
    std::cout << "\tUse sim input  : \t" << bUseSim << std::endl;
    std::cout << "\tSim input file : \t" << siminput << std::endl;
    if (argc > 9)
        std::cout << "\tRun over events from " << evStart << " to " << nEvents << " (only if input is used)" << std::endl;

    detector det = kJFoCal;

    TFile *fOut = new TFile(outFileName, "RECREATE");

    Pythia pythia;

    AliJHMRHist *fHistos = new AliJHMRHist();
    fHistos->CreateHistos(fOut, det, bUseLeading);

    AliJHMRPythiaCatalyst *fCatalyst = new AliJHMRPythiaCatalyst(pythia.event, fHistos);
    AliJHMRGeantCatalyst *fCatalystG = new AliJHMRGeantCatalyst(siminput, fHistos);

    AliJHMRCorr *fCorr = new AliJHMRCorr(fHistos, det, bUseLeading, bUseSim);

    TClonesArray *arrPhotonFor      = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Matched     = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Real        = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Peak        = new TClonesArray("AliJBaseTrack", 1500);
    TClonesArray *arrPi0Side        = new TClonesArray("AliJBaseTrack", 1500);

    TClonesArray* arrPi0PeakMixed[poolsize];
    TClonesArray* arrPi0SideMixed[poolsize];

    for (int ipool = 0; ipool < poolsize; ipool++) {
        arrPi0PeakMixed[ipool] = new TClonesArray("AliJBaseTrack", 1500);
        arrPi0SideMixed[ipool] = new TClonesArray("AliJBaseTrack", 1500);
    }

    if (bUseSim) {
        if (!fCatalystG->LoadInput()) return 1;
        if (nEvents==-1)
            nEvents = fCatalystG->GetNumberOfEvents();
    } else {
        // Initialise pythia
        pythia.readFile(pythiaSettings.Data());
        pythia.readString("Random:setSeed = on");
        pythia.readString(Form("Random:seed = %d", seed));
        pythia.init();
        evStart = 0;
        nEvents = pythia.mode("Main:numberOfEvents");
    }

    std::cout << "Number of events : " << nEvents - evStart << std::endl;

    //
    // Loop over events
    //
    for ( int iEvent = evStart; iEvent < nEvents; ++iEvent ) {

        fHistos->hCounter->Fill(0.5); // number of events

        if (bDebugOn)
            std::cout << "\nEvent " << iEvent << std::endl;

        if (bUseSim) { // Get stuff from detector simulation file
            if (iEvent%10000==0) cout << "event " << iEvent << "/" << nEvents << endl;
            fCatalystG->GetEvent(iEvent);
            fCatalystG->InitializeEvent();
            fCatalystG->GetParticles();
            fCatalystG->GetClusters();
            arrPi0Real   = fCatalystG->GetParticleList(kJPi0);
            arrPhotonFor = fCatalystG->GetParticleList(kJCluster);
        } else { // Get stuff from pythia
            if ( !pythia.next() ) continue;
            fCatalyst->InitializeEvent();
            fCatalyst->GetParticles(det);
            fHistos->hTotCrossSec->Fill(pythia.info.sigmaGen());
            arrPhotonFor = fCatalyst->GetParticleList(kJDecayPhoton);
            arrPi0Real   = fCatalyst->GetParticleList(kJPi0);
            fCorr->SmearEnergies(arrPhotonFor);
            fCorr->SetX1X2(pythia.info.x1(), pythia.info.x2());
        }
        int nTrueFromPeak = fCorr->ReconstructPions(arrPhotonFor, arrPi0Peak, det, 1);
        fCorr->ReconstructPions(arrPhotonFor, arrPi0Side, det, 0);

        // Do matching for clusters
        if (bUseSim) {
            fCatalystG->GetPtMatchedClusters(arrPi0Peak);
            arrPi0Matched = fCatalystG->GetParticleList(kJPtMatchedPi0);
            fHistos->FillMathingInformation(arrPi0Peak, arrPi0Matched);
        }

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
        fCorr->GetTriggAssocLists(arrPi0Real, listTriggReal, listAssocReal, binsWithTriggReal, 0);
        fCorr->GetTriggAssocLists(arrPi0Peak, listTriggPeak, listAssocPeak, binsWithTriggPeak, 0);
        //fCorr->GetTriggAssocLists(arrPi0Matched, listTriggPeak, listAssocPeak, binsWithTriggPeak, 0);
        fCorr->GetTriggAssocLists(arrPi0Side, listTriggSide, listAssocSide, binsWithTriggSide, 0);

        if (bUseLeading) {
            int isPeakTriggLarger = fCorr->GetLargerTrigg(arrPi0Peak, listTriggPeak, arrPi0Side, listTriggSide);
            if (isPeakTriggLarger==1)
                fCorr->FillPionMassesLeading(arrPhotonFor, arrPi0Peak, listTriggPeak, det, isPeakTriggLarger);
            if (isPeakTriggLarger==0)
                fCorr->FillPionMassesLeading(arrPhotonFor, arrPi0Side, listTriggSide, det, isPeakTriggLarger);
        } else {
            fCorr->FillPionMasses(arrPhotonFor, binsWithTriggPeak, binsWithTriggSide, det);
        }

        //fCorr->FillPionMassesTrue(arrPi0Real, binsWithTriggReal, det);
        fCorr->FillRealTriggers(arrPi0Real, listTriggReal);
        fCorr->FillAsymmetry(arrPhotonFor, det);

        fCorr->DoCorrelations(arrPi0Real, arrPhotonFor, listTriggReal, listAssocReal, fHistos->hCorrFor, 1, 0, 0); // last three values: bTrueCorr, bMassWindow, bUseWeight
        fCorr->DoCorrelations(arrPi0Peak, arrPhotonFor, listTriggPeak, listAssocPeak, fHistos->hCorrMassMass, 0, 1, 1);

        if (bUseLeading) {
            int isPeakTriggLarger = fCorr->GetLargerTrigg(arrPi0Peak, listTriggPeak, arrPi0Side, listTriggSide);
            if (isPeakTriggLarger) {
                fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, fHistos->hCorrMassSide, 0);
            } else {
                fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, fHistos->hCorrSideMass, 0);
                fCorr->DoCorrelations(arrPi0Side, arrPhotonFor, listTriggSide, listAssocSide, fHistos->hCorrSideSide, 0, 0, 0);
            }
        } else {
            fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0Side, listAssocSide, fHistos->hCorrMassSide, 1);
            fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0Peak, listAssocPeak, fHistos->hCorrSideMass, 1);
            fCorr->DoCorrelations(arrPi0Side, arrPhotonFor, listTriggSide, listAssocSide, fHistos->hCorrSideSide, 0, 0, 1);
        }

        // Construct & save true correlation components f_SS, f_SB, f_BS, f_BB
        fCorr->ConstructTrueCorrComponents(arrPi0Peak, arrPhotonFor, listTriggPeak, listAssocPeak, 0);
        fCorr->ClearPhotonPairVector();
        fCorr->ClearSidebandPairVector();

        // Mixed event : take triggers from this event, associated from previous
        if (poolsize > 0) {
            if (iEvent >= poolsize) {
                for (int ipool = 0; ipool < poolsize; ipool++) {
                    std::vector<int> listTriggPeakMixed, listTriggSideMixed, listAssocPeakMixed, listAssocSideMixed;
                    int binsWithTriggPeakMixed[NTRIGGBINS+1] = {0}, binsWithTriggSideMixed[NTRIGGBINS+1] = {0};
                    fCorr->GetTriggAssocLists(arrPi0PeakMixed[ipool], listTriggPeakMixed, listAssocPeakMixed, binsWithTriggPeakMixed, 0);
                    fCorr->GetTriggAssocLists(arrPi0SideMixed[ipool], listTriggSideMixed, listAssocSideMixed, binsWithTriggSideMixed, 0);
                    fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0PeakMixed[ipool], listAssocPeakMixed, fHistos->hCorrMassMassMixed, 0);
                    fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0SideMixed[ipool], listAssocSideMixed, fHistos->hCorrSideSideMixed, 0);
                    fCorr->DoCorrelations(arrPi0Peak, listTriggPeak, arrPi0SideMixed[ipool], listAssocSideMixed, fHistos->hCorrMassSideMixed, 0);
                    fCorr->DoCorrelations(arrPi0Side, listTriggSide, arrPi0PeakMixed[ipool], listAssocPeakMixed, fHistos->hCorrSideMassMixed, 0);

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
