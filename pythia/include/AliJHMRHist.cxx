#include "AliJHMRHist.h"

void AliJHMRHist::CreateHistos(TFile *output, detector labelDet) {

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

    // Correlation and mass histograms
    dirMasses = output->mkdir("Masses");
    dirCorrMid = output->mkdir("CorrMid");
    dirCorrFor = output->mkdir("CorrFor");
    dirCorrChargedMid = output->mkdir("CorrChargedMid");
    dirCorrChargedFor = output->mkdir("CorrChargedFor");
    dirCorrMassMass = output->mkdir("CorrMassMass");
    dirCorrMassSide = output->mkdir("CorrMassSide");
    dirCorrSideMass = output->mkdir("CorrSideMass");
    dirCorrSideSide = output->mkdir("CorrSideSide");

    for (int i = 0; i < NTRIGGBINS; i++) {
        double tlow = triggPt[i];
        double tupp = triggPt[i+1];
        dirMasses->cd();
        hPi0MassTrigg[i] = new TH1D(Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), 360, 0.0, 720.0);
        hPi0MassTrigg[i]->Sumw2();

        for (int j = 0; j < NASSOCBINS; j++) {
            double alow = assocPt[j];
            double aupp = assocPt[j+1];

            if (tlow < aupp) continue;

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

        }
    }
}

void AliJHMRHist::FillPtEta(particleType itype, TClonesArray * arrParticles)
{
    int nent = arrParticles->GetEntriesFast();
    for (int ient = 0; ient < nent; ient++) {
        TLorentzVector *lv = (TLorentzVector*)arrParticles->At(ient);
        double pT = lv->Pt();
        double eta = lv->Eta();
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
            hRecPionPt->Fill(pT);
            hRecPionEta->Fill(eta);
        }
    }
}
