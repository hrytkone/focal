#include "Efficiency.h"

void Efficiency(TString inputfile)
{
    fGeom = AliFOCALGeometry::GetInstance("../geometry.txt");
    cout << "FoCal front at " << (double)fGeom->GetFOCALZ0() - (double)fGeom->GetFOCALSizeZ()/2 << endl;

    if (!LoadInput(inputfile)) {
        cout << "File not found! " << inputfile << endl;
        return;
    }

    InitOutput();

    int nev = fTree->GetEntries();
    cout << "Processing " << nev << " events" << endl;
    for (int iev=0; iev<nev; iev++) {

        if (iev%1000==0) cout << "event " << iev << "/" << nev << endl;
        //cout << "\n================================" << endl;
        //cout << "Event " << iev << endl;

        fTree->GetEntry(iev);

        int nclust = fClusters->GetEntries();
        for (int iclust=0; iclust<nclust; iclust++) {
            AliJHMRCluster *clust = (AliJHMRCluster*)fClusters->At(iclust);
            double clustEta = clust->GetEta();
            double clustPhi = clust->GetPhi();
            double clustX = clust->GetX();
            double clustY = clust->GetY();
            double clustZ = clust->GetZ();
            double clustE = clust->GetE();
            double clustEHCAL = clust->GetEHCAL();
            double clustPt = clust->GetPt();

            hEPhotonCluster->Fill(clustE);
            hPhiEtaGamma->Fill(clustPhi, clustEta);
            hPhiThetaGamma->Fill(clustPhi, 2.*TMath::ATan(TMath::Exp(-clustEta)));
            hXYGamma->Fill(clustX, clustY);

            lvParticle.SetPtEtaPhiM(clustPt, clustEta, clustPhi, 0.);
            AliJBaseTrack track( lvParticle );
            new((*fTrackClusters)[fTrackClusters->GetEntriesFast()]) AliJBaseTrack(track);
        }

        FillTruePions();
        FillRecPions(fTrackClusters);
        fTrackClusters->Clear("C");

    }
    cout << "NTRUE = " << ntrue << endl;
    fOut->cd();
    fOut->Write("", TObject::kOverwrite);
}

Int_t LoadInput(TString inputfile)
{
    fIn = TFile::Open(inputfile.Data(), "READ");
    if (!fIn) return 0;

    fTree = (TTree*)fIn->Get("Run");
    fTree->SetBranchAddress("Events", &fEvent);
    fBrtrack = fTree->GetBranch("Events.fTracks");
    fBrtrack->SetAddress(&fTracks);
    fBrcluster = fTree->GetBranch("Events.fClusters");
    fBrcluster->SetAddress(&fClusters);
    return 1;
}

void InitOutput()
{
    fOut = TFile::Open(Form("%s.root", outputname.Data()), "RECREATE");
    cout << "\npt=[ ";
    for (int i = 0; i <= nPtBin; i++) {
        pt[i] = limMin*exp(i*logBW);
        cout << pt[i] << " ";
    }
    cout << "]\neta=[ ";
    for (int i = 0; i <= nEtaBin; i++) {
        eta[i] = etamin + i*etaBW;
        cout << eta[i] << " ";
    }
    cout << "]\n" << endl;

    hEtaPtTrue = new TH2D("hEtaPtTrue", "hEtaPtTrue", nEtaBin, eta, nPtBin, pt); hEtaPtTrue->Sumw2();
    hEtaPtRec = new TH2D("hEtaPtRec", "hEtaPtRec", nEtaBin, eta, nPtBin, pt); hEtaPtRec->Sumw2();
    hEtaPtRec_match = new TH2D("hEtaPtRec_match", "hEtaPtRec_match", nEtaBin, eta, nPtBin, pt); hEtaPtRec_match->Sumw2();
    hEtaPtRec_match_truept = new TH2D("hEtaPtRec_match_truept", "hEtaPtRec_match_truept", nEtaBin, eta, nPtBin, pt); hEtaPtRec_match_truept->Sumw2();
    hEtaETrue = new TH2D("hEtaETrue", "hEtaETrue", nEtaBin, eta, 200, 0., 200.);
    hEtaERec = new TH2D("hEtaERec", "hEtaERec", nEtaBin, eta, 200, 0., 200.);
    hEtaERec_match = new TH2D("hEtaERec_match", "hEtaERec_match", nEtaBin, eta, 200, 0., 200.);
    hEPhotonTrue = new TH1D("hEPhotonTrue", "hEPhotonTrue", 200, 0., 200.);
    hEPhotonCluster = new TH1D("hEPhotonCluster", "hEPhotonCluster", 200, 0., 200.);

    hPtMass = new TH2D("hPtMass", "hPtMass", nPtBin, limMin, limMax, 350, 0., 700.);
    hPtMass_match = new TH2D("hPtMass_match", "hPtMass_match", nPtBin, limMin, limMax, 350, 0., 700.);
    hEtaMass = new TH2D("hEtaMass", "hEtaMass", nEtaBin, eta, 350, 0., 700.);
    hEtaMass_match = new TH2D("hEtaMass_match", "hEtaMass_match", nEtaBin, eta, 350, 0., 700.);

    hPhiEtaTrue = new TH2D("hPhiEtaTrue", "hPhiEtaTrue", nPhiBin, phimin, phimax, nEtaBin, eta);
    hPhiEta = new TH2D("hPhiEta", "hPhiEta", nPhiBin, phimin, phimax, nEtaBin, eta);
    hPhiEta_match = new TH2D("hPhiEta_match", "hPhiEta_match", nPhiBin, phimin, phimax, nEtaBin, eta);
    hPhiTheta = new TH2D("hPhiTheta", "hPhiTheta", nPhiBin, phimin, phimax, nThetaBin, thetamin, thetamax);
    hPhiTheta_match  = new TH2D("hPhiTheta_match", "hPhiTheta_match", nPhiBin, phimin, phimax, nThetaBin, thetamin, thetamax);
    hXY = new TH2D("hXY", "hXY", nXYBin, xymin, xymax, nXYBin, xymin, xymax);
    hXY_match = new TH2D("hXY_match", "hXY_match", nXYBin, xymin, xymax, nXYBin, xymin, xymax);

    hPhiEtaGamma = new TH2D("hPhiEtaGamma", "hPhiEtaGamma", nPhiBin, phimin, phimax, nEtaBin, eta);
    hPhiThetaGamma = new TH2D("hPhiThetaGamma", "hPhiThetaGamma", nPhiBin, phimin, phimax, nThetaBin, thetamin, thetamax);
    hXYGamma = new TH2D("hXYGamma", "hXYGamma", nXYBin, xymin, xymax, nXYBin, xymin, xymax);
}

void FillTruePions()
{
    int ntrack = fTracks->GetEntries();
    for (int itrack=0; itrack<ntrack; itrack++) {
        AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(itrack);
        double trEta = tr->Eta();
        double trPhi = tr->Phi();
        double trPt = tr->Pt();
        double trE = tr->E();
        double trPid = tr->GetID();
        double trMomType = tr->GetMotherType();
        if (trEta < etamin || trEta > etamax) continue;
        if (trPt < limMin || trPt > limMax) continue;
        if (trPid==111) {
            hEtaPtTrue->Fill(trEta, trPt);
            hEtaETrue->Fill(trEta, trE);
            hPhiEtaTrue->Fill(trPhi, trEta);
            ntrue++;
        }
        if (trPid==22 && trMomType==111 && tr->IsPrimary())
            hEPhotonTrue->Fill(trE);
    }
}

void FillRecPions(TClonesArray *clusters)
{
    int nclust = clusters->GetEntriesFast();
    if (nclust < 2) return;
    double z = 700.;

    double closestPt = -1.;
    double closestEta = -1.;
    double closestTheta = -1.;
    double closestPhi = -1.;
    double closestE = -1.;
    double closestMass = -1.;
    double closestAsym = -1.;
    double massdiff = 10000.;
    for (int i = 1; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)clusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)clusters->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(lv1, lv2);
            if (lvSum.Pt()<limMin || lvSum.Pt()>limMax) continue;
            if (lvSum.Eta()<etamin || lvSum.Eta()>etamax) continue;
            double mass = 1000.*lvSum.M();
            if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())<asymcut) {
                hPtMass->Fill(lvSum.Pt(), mass);
                hEtaMass->Fill(lvSum.Eta(), mass);
                if (mass > 100. && mass < 150.) {
                    hPhiEta->Fill(lvSum.Phi(), lvSum.Eta());
                    hPhiTheta->Fill(lvSum.Phi(), lvSum.Theta());
                    hXY->Fill(z*TMath::Tan(lvSum.Theta())*TMath::Cos(lvSum.Phi()), z*TMath::Tan(lvSum.Theta())*TMath::Sin(lvSum.Phi()));
                    hEtaERec->Fill(lvSum.Eta(), lvSum.E());
                    hEtaPtRec->Fill(lvSum.Eta(), lvSum.Pt());
                }
            }

            if (TMath::Abs(mass-135.0) < massdiff) {
                massdiff = TMath::Abs(mass-135.0);
                closestMass = mass;
                closestPt = lvSum.Pt();
                closestEta = lvSum.Eta();
                closestPhi = lvSum.Phi();
                closestE = lvSum.E();
                closestTheta = lvSum.Theta();
                closestAsym = TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E());
            }
        }
    }
    AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(0);
    double truePt = tr->Pt();
    
    if (closestMass > 100. && closestMass < 150. && closestAsym < asymcut) {
        if (closestEta>etamin && closestEta<etamax)
            hPtMass_match->Fill(closestPt, closestMass);
        if (closestPt>limMin && closestPt<limMax)
            hEtaMass_match->Fill(closestEta, closestMass); 
        hPhiEta_match->Fill(closestPhi, closestEta);
        hPhiTheta_match->Fill(closestPhi, closestTheta);
        hXY_match->Fill(z*TMath::Tan(closestTheta)*TMath::Cos(closestPhi), z*TMath::Tan(closestTheta)*TMath::Sin(closestPhi));
        hEtaERec_match->Fill(closestEta, closestE);
        hEtaPtRec_match->Fill(closestEta, closestPt);
        hEtaPtRec_match_truept->Fill(closestEta, truePt);
    }

}

int GetBin(double arr[], int nArr, double val)
{
    for (int i=0; i<nArr; i++) {
        if (arr[i]<=val && val<arr[i+1]) return i;
    }
    return -1;
}

AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    AliJBaseTrack lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}
