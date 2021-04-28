#include "AliJHMRHist.h"

void CreateHistos() {
    // Basic histograms
    TH1D *hCounter = new TH1D("hCounter", "hCounter", 10, 0, 10);
    TH1D *hRealTriggCounter = new TH1D("hRealTriggCounter", "hRealtriggCounter", 10, 0, 10);
    
    for (int i = 0; i <= nIncPtBin; i++) logBinsX[i] = limMin*exp(i*logBW);
    
    TH1D *hPionPt = new TH1D("hPionPt", "hPionPt", nIncPtBin, logBinsX); hPionPt->Sumw2();
    TH1D *hPionPtFor = new TH1D("hPionPtFor", "hPionPtFor", nIncPtBin, logBinsX); hPionPtFor->Sumw2();
    TH1D *hPionPtForDetected = new TH1D("hPionPtForDetected", "hPionPtForDetected", nIncPtBin, logBinsX); hPionPtForDetected->Sumw2();
    TH1D *hPionPtMid = new TH1D("hPionPtMid", "hPionPtMid", nIncPtBin, logBinsX); hPionPtMid->Sumw2();
    
    TH1D *hChargedHadronPt = new TH1D("hChargedHadronPt", "hChargedHadronPt", nIncPtBin, logBinsX); hChargedHadronPt->Sumw2();
    TH1D *hChargedHadronPtFor = new TH1D("hChargedHadronPtFor", "hChargedHadronPtFor", nIncPtBin, logBinsX); hChargedHadronPtFor->Sumw2();
    TH1D *hChargedHadronPtMid = new TH1D("hChargedHadronPtMid", "hChargedHadronPtMid", nIncPtBin, logBinsX); hChargedHadronPtMid->Sumw2();

    TH1D *hPhotonPt = new TH1D("hPhotonPt", "hPhotonPt", nIncPtBin, logBinsX); hPhotonPt->Sumw2();
    TH1D *hPhotonPtFor = new TH1D("hPhotonPtFor", "hPhotonPtFor", nIncPtBin, logBinsX); hPhotonPtFor->Sumw2();
    TH1D *hPhotonPtMid = new TH1D("hPhotonPtMid", "hPhotonPtMid", nIncPtBin, logBinsX); hPhotonPtMid->Sumw2();
    TH1D *hPhotonEnergyReal = new TH1D("hPhotonEnergyReal", "hPhotonEnergyReal", nPhotonEnergyBin, limPhotonEnergyMin, limPhotonEnergyMax); hPhotonEnergyReal->Sumw2();
    TH1D *hPhotonEnergy = new TH1D("hPhotonEnergy", "hPhotonEnergy", nPhotonEnergyBin, limPhotonEnergyMin, limPhotonEnergyMax); hPhotonEnergy->Sumw2();
    
    TH1D *hPionEta = new TH1D("hPionEta", "hPionEta", nIncEtaBin, -incEtaRange/2., incEtaRange/2.); hPionEta->Sumw2();
    TH1D *hChargedHadronEta = new TH1D("hChargedHadronEta", "hChargedHadronEta", nIncEtaBin, -incEtaRange/2., incEtaRange/2.); hChargedHadronEta->Sumw2();

    // Correlation and mass histograms
    TDirectory *dirMasses = fOut->mkdir("Masses");
    TDirectory *dirCorrMid = fOut->mkdir("CorrMid");
    TDirectory *dirCorrFor = fOut->mkdir("CorrFor");
    TDirectory *dirCorrChargedMid = fOut->mkdir("CorrChargedMid");
    TDirectory *dirCorrChargedFor = fOut->mkdir("CorrChargedFor");
    TDirectory *dirCorrMassMass = fOut->mkdir("CorrMassMass");
    TDirectory *dirCorrMassSide = fOut->mkdir("CorrMassSide");
    TDirectory *dirCorrSideMass = fOut->mkdir("CorrSideMass");
    TDirectory *dirCorrSideSide = fOut->mkdir("CorrSideSide");

    TH1D *hPi0MassTrigg[nTriggBins];
    TH1D *hPi0MassAssocPeak[nTriggBins][nAssocBins];
    TH1D *hPi0MassAssocSide[nTriggBins][nAssocBins];
 
    TH2D *hCorrMid[nTriggBins][nAssocBins];
    TH2D *hCorrFor[nTriggBins][nAssocBins];
    TH2D *hCorrChargedMid[nTriggBins][nAssocBins];
    TH2D *hCorrChargedFor[nTriggBins][nAssocBins];
    
    TH2D *hCorrMassMass[nTriggBins][nAssocBins];
    TH2D *hCorrMassSide[nTriggBins][nAssocBins];
    TH2D *hCorrSideMass[nTriggBins][nAssocBins];
    TH2D *hCorrSideSide[nTriggBins][nAssocBins];

    for (int i = 0; i < nTriggBins; i++) {
        double tlow = triggPt[i];
        double tupp = triggPt[i+1];
        dirMasses->cd();
        hPi0MassTrigg[i] = new TH1D(Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), Form("hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp), 360, 0.0, 720.0);
        hPi0MassTrigg[i]->Sumw2();

        for (int j = 0; j < nAssocBins; j++) {
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
            hCorrMid[i][j] = new TH2D(Form("hCorrMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -2.*etaTrackerRange, 2.*etaTrackerRange);
            hCorrMid[i][j]->Sumw2();
            
            dirCorrFor->cd();
            hCorrFor[i][j] = new TH2D(Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp),
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrFor[i][j]->Sumw2();
            
            dirCorrChargedMid->cd();
            hCorrChargedMid[i][j] = new TH2D(Form("hCorrChargedMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrChargedMid[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinTracker, -2.*etaTrackerRange, 2.*etaTrackerRange);
            hCorrChargedMid[i][j]->Sumw2();
            
            dirCorrChargedFor->cd();
            hCorrChargedFor[i][j] = new TH2D(Form("hCorrChargedFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrChargedFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrChargedFor[i][j]->Sumw2();
            
            dirCorrMassMass->cd();
            hCorrMassMass[i][j] = new TH2D(Form("hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrMassMass[i][j]->Sumw2();

            dirCorrMassSide->cd();
            hCorrMassSide[i][j] = new TH2D(Form("hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrMassSide[i][j]->Sumw2();

            dirCorrSideMass->cd();
            hCorrSideMass[i][j] = new TH2D(Form("hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrSideMass[i][j]->Sumw2();

            dirCorrSideSide->cd();
            hCorrSideSide[i][j] = new TH2D(Form("hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), Form("hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), 
                                    nPhiBin, deltaPhiMin, deltaPhiMax, nEtaBinFocal, -etaFocalRange, etaFocalRange);
            hCorrSideSide[i][j]->Sumw2();
            
        }
    }
}

