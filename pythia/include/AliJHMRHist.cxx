#include "AliJHMRHist.h"

void AliJHMRHist::CreateHistos(TFile *output, detector labelDet, bool bUseLeading) {

    double etaRange = 0;
    int nEtaBin = 0;
    if (labelDet==kJFoCal) {
        etaRange = etaFocalRange;
        nEtaBin = nEtaBinFocal;
    }
    if (labelDet==kJSTAR) {
        etaRange = etaSTARRange;
        nEtaBin = nEtaBinSTAR;
    }

    // Basic histograms
    hCounter = new TH1D("hCounter", "hCounter", 10, 0, 10);
    hRealTriggCounter = new TH1D("hRealTriggCounter", "hRealtriggCounter", 10, 0, 10);

    for (int i = 0; i <= NINCPTBIN; i++) logBinsX[i] = limMin*exp(i*logBW);

    hPionPt = new TH1D("hPionPt", "hPionPt", NINCPTBIN, logBinsX); hPionPt->Sumw2();
    hRecPionPt = new TH1D("hRecPionPt", "hRecPionPt", NINCPTBIN, logBinsX); hRecPionPt->Sumw2();
    hChargedHadronPt = new TH1D("hChargedHadronPt", "hChargedHadronPt", NINCPTBIN, logBinsX); hChargedHadronPt->Sumw2();
    hPhotonPt = new TH1D("hPhotonPt", "hPhotonPt", NINCPTBIN, logBinsX); hPhotonPt->Sumw2();

    hPionEta = new TH1D("hPionEta", "hPionEta", NINCETABIN, -incEtaRange/2., incEtaRange/2.); hPionEta->Sumw2();
    hRecPionEta = new TH1D("hPionEta", "hPionEta", NINCETABIN, -incEtaRange/2., incEtaRange/2.); hRecPionEta->Sumw2();
    hPhotonEta = new TH1D("hPhotonEta", "hPhotonEta", NINCETABIN, -incEtaRange/2., incEtaRange/2.); hPhotonEta->Sumw2();
    hChargedHadronEta = new TH1D("hChargedHadronEta", "hChargedHadronEta", NINCETABIN, -incEtaRange/2., incEtaRange/2.); hChargedHadronEta->Sumw2();

    hEnergyAsymTrue = new TH1D("hEnergyAsymTrue", "hEnergyAsymTrue", 100, 0., 1.); hEnergyAsymTrue->Sumw2();
    hEnergyAsymRec = new TH1D("hEnergyAsymRec", "hEnergyAsymRec", 100, 0., 1.); hEnergyAsymRec->Sumw2();

    // Correlation and mass histograms
    dirMasses = output->mkdir("Masses");
    dirMassComponents = output->mkdir("MassComponents");
    dirCorrMid = output->mkdir("CorrMid");
    dirCorrFor = output->mkdir("CorrFor");
    dirCorrMeas = output->mkdir("CorrMeas");
    dirCorrReject = output->mkdir("CorrReject");
    dirCorrChargedMid = output->mkdir("CorrChargedMid");
    dirCorrChargedFor = output->mkdir("CorrChargedFor");
    dirCorrMassMass = output->mkdir("CorrMassMass");
    dirCorrMassSide = output->mkdir("CorrMassSide");
    dirCorrSideMass = output->mkdir("CorrSideMass");
    dirCorrSideSide = output->mkdir("CorrSideSide");
    dirTrueComponents = output->mkdir("TrueComponents");
    dirBackgroundSources = output->mkdir("BackgroundSources");
    dirBackgroundEnergies = output->mkdir("BackgroundEnergies");

    dirMassComponents->cd();
    hMassTrue = new TH1D("hMassTrue", "Invariant mass, photons from same mother", 360, 0.0, 720.0);
    hMassTrue->Sumw2();
    hMassDecay = new TH1D("hMassDecay", "Invariant mass, photons from pi0 decay but from different mothers", 360, 0.0, 720.0);
    hMassDecay->Sumw2();
    hMassMix = new TH1D("hMassMix", "Invariant mass, only one photon is from pi0 decay", 360, 0.0, 720.0);
    hMassMix->Sumw2();
    hMassNotDecay = new TH1D("hMassNotDecay", "Invariant mass, either from the photons is from pi0 decay", 360, 0.0, 720.0);
    hMassNotDecay->Sumw2();

    for (int i = 0; i < NTRIGGBINS; i++) {
        double tlow = triggPt[i];
        double tupp = triggPt[i+1];
        dirMasses->cd();
        hPi0MassTrigg[i] = new TH1D(Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), 360, 0.0, 720.0);
        hPi0MassTrigg[i]->Sumw2();

        dirBackgroundEnergies->cd();
        hEnergyMassBgTrigg[i] = new TH1D(Form("hEnergyMassBgTrigg[%4.1f,%4.1f]",tlow,tupp), Form("hEnergyMassBgTrigg[%4.1f,%4.1f]",tlow,tupp), 200, 0., 100.);
        hEnergyMassBgTrigg[i]->Sumw2();
        hEnergySidebandTrigg[i] = new TH1D(Form("hEnergySidebandTrigg[%4.1f,%4.1f]",tlow,tupp), Form("hEnergySidebandTrigg[%4.1f,%4.1f]",tlow,tupp), 200, 0., 100.);
        hEnergySidebandTrigg[i]->Sumw2();

        for (int j = 0; j < NASSOCBINS; j++) {
            double alow = assocPt[j];
            double aupp = assocPt[j+1];

            if (tlow < aupp && !bUseLeading) continue;

            dirMasses->cd();
            hPi0MassAssocPeak[i][j] = new TH1D(Form("hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                               Form("hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 360, 0.0, 720.0);
            hPi0MassAssocPeak[i][j]->Sumw2();
            hPi0MassAssocSide[i][j] = new TH1D(Form("hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                               Form("hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 360, 0.0, 720.0);
            hPi0MassAssocSide[i][j]->Sumw2();

            dirCorrMid->cd();
            hCorrMid[i][j] = new TH2D(Form("hCorrMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      Form("hCorrMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -2.*etaTrackerRange, 2.*etaTrackerRange);
            hCorrMid[i][j]->Sumw2();

            dirCorrFor->cd();
            hCorrFor[i][j] = new TH2D(Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrFor[i][j]->Sumw2();

            dirCorrMeas->cd();
            hCorrMeas[i][j] = new TH2D(Form("hCorrMeas[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      Form("hCorrMeas[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrMeas[i][j]->Sumw2();

            dirCorrReject->cd();
            hCorrRejectMassMass[i][j] = new TH2D(Form("hCorrRejectMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      Form("hCorrRejectMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrRejectMassMass[i][j]->Sumw2();
            hCorrRejectSideSide[i][j] = new TH2D(Form("hCorrRejectSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      Form("hCorrRejectSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                      nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrRejectSideSide[i][j]->Sumw2();


            dirCorrChargedMid->cd();
            hCorrChargedMid[i][j] = new TH2D(Form("hCorrChargedMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                             Form("hCorrChargedMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                             nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -2.*etaTrackerRange, 2.*etaTrackerRange);
            hCorrChargedMid[i][j]->Sumw2();

            dirCorrChargedFor->cd();
            hCorrChargedFor[i][j] = new TH2D(Form("hCorrChargedFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                             Form("hCorrChargedFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                             nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrChargedFor[i][j]->Sumw2();

            dirCorrMassMass->cd();
            hCorrMassMass[i][j] = new TH2D(Form("hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrMassMass[i][j]->Sumw2();
            hCorrMassMassMixed[i][j] = new TH2D(Form("hCorrMassMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                Form("hCorrMassMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrMassMassMixed[i][j]->Sumw2();

            dirCorrMassSide->cd();
            hCorrMassSide[i][j] = new TH2D(Form("hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrMassSide[i][j]->Sumw2();
            hCorrMassSideMixed[i][j] = new TH2D(Form("hCorrMassSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                Form("hCorrMassSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrMassSideMixed[i][j]->Sumw2();

            dirCorrSideMass->cd();
            hCorrSideMass[i][j] = new TH2D(Form("hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSideMass[i][j]->Sumw2();
            hCorrSideMassMixed[i][j] = new TH2D(Form("hCorrSideMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                Form("hCorrSideMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSideMassMixed[i][j]->Sumw2();

            dirCorrSideSide->cd();
            hCorrSideSide[i][j] = new TH2D(Form("hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSideSide[i][j]->Sumw2();
            hCorrSideSideMixed[i][j] = new TH2D(Form("hCorrSideSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                Form("hCorrSideSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                                nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSideSideMixed[i][j]->Sumw2();

            dirTrueComponents->cd();
            hCorrSignalSignal[i][j] = new TH2D(Form("hCorrSignalSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrSignalSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSignalSignal[i][j]->Sumw2();
            hCorrSignalBg[i][j] = new TH2D(Form("hCorrSignalBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrSignalBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSignalBg[i][j]->Sumw2();
            hCorrBgSignal[i][j] = new TH2D(Form("hCorrBgSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrBgSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrBgSignal[i][j]->Sumw2();
            hCorrBgBg[i][j] = new TH2D(Form("hCorrBgBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrBgBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrBgBg[i][j]->Sumw2();


            dirBackgroundSources->cd();
            hCorrBgBgDecay[i][j] = new TH2D(Form("hCorrBgBgDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrBgBgDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrBgBgDecay[i][j]->Sumw2();
            hCorrBgBgMix[i][j] = new TH2D(Form("hCorrBgBgMix[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrBgBgMix[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrBgBgMix[i][j]->Sumw2();
            hCorrBgBgNotDecay[i][j] = new TH2D(Form("hCorrBgBgNotDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrBgBgNotDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrBgBgNotDecay[i][j]->Sumw2();

            hCorrSideSideDecay[i][j] = new TH2D(Form("hCorrSideSideDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrSideSideDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSideSideDecay[i][j]->Sumw2();
            hCorrSideSideMix[i][j] = new TH2D(Form("hCorrSideSideMix[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrSideSideMix[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSideSideMix[i][j]->Sumw2();
            hCorrSideSideNotDecay[i][j] = new TH2D(Form("hCorrSideSideNotDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hCorrSideSideNotDecay[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBin, -etaRange, etaRange);
            hCorrSideSideNotDecay[i][j]->Sumw2();

            dirBackgroundEnergies->cd();
            hEnergyMassBgAssoc[i][j] = new TH1D(Form("hEnergyMassBgAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           Form("hEnergyMassBgAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                           200, 0., 1000.);
            hEnergyMassBgAssoc[i][j]->Sumw2();
            hEnergySidebandAssoc[i][j] = new TH1D(Form("hEnergySidebandAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                             Form("hEnergySidebandAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                             200, 0., 1000.);
            hEnergySidebandAssoc[i][j]->Sumw2();
        }
    }
}

void AliJHMRHist::FillPtEta(particleType itype, TClonesArray * arrParticles)
{
    int nent = arrParticles->GetEntriesFast();
    for (int ient = 0; ient < nent; ient++) {
        AliJBaseTrack *lv = (AliJBaseTrack*)arrParticles->At(ient);
        double pT = lv->Pt();
        double eta = lv->Eta();
        int isPion = lv->GetLabel();
        if (itype==kJHadron) {
            hChargedHadronPt->Fill(pT);
            hChargedHadronEta->Fill(eta);
        }

        if (itype==kJPi0) {
            hPionPt->Fill(pT);
            hPionEta->Fill(eta);
        }

        if (itype==kJDecayPhoton) {
            hPhotonPt->Fill(pT);
            hPhotonEta->Fill(eta);
        }

        if (itype==kJRecPi0) {
            if (isPion) {
                hRecPionPt->Fill(pT);
                hRecPionEta->Fill(eta);
            }
        }
    }
}
