#include "Efficiency.h"

void Efficiency(TString inputfile)
{
    fGeom = AliFOCALGeometry::GetInstance("../geometry.txt");

    if (!LoadInput(inputfile)) {
        cout << "File not found! " << inputfile << endl;
        return;
    }

    InitOutput();

    //int nev = 1000;
    int nev = fTree->GetEntries();
    cout << "Processing " << nev << " events" << endl;
    for (int iev=0; iev<nev; iev++) {

        if (iev%100==0) cout << "event " << iev << "/" << nev << endl;
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
            double clustE = clust->GetE();
            double clustEHCAL = clust->GetEHCAL();
            double clustPt = clust->GetPt();

            hEPhotonCluster->Fill(clustE);

            lvParticle.SetPtEtaPhiM(clustPt, clustEta, clustPhi, 0.);
            AliJBaseTrack track( lvParticle );
            new((*fTrackClusters)[fTrackClusters->GetEntriesFast()]) AliJBaseTrack(track);
        }

        FillTruePions();
        FillRecPions();
        fTrackClusters->Clear("C");

    }
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
    for (int i = 0; i <= nIncPtBin; i++) logBinsX[i] = limMin*exp(i*logBW);
    hEtaPtTrue = new TH2D("hEtaPtTrue", "hEtaPtTrue", nEtaBin, etamin, etamax, nIncPtBin, logBinsX);
    hEtaETrue = new TH2D("hEtaETrue", "hEtaETrue", nEtaBin, etamin, etamax, 200, 0., 200.);
    hEtaPtRec = new TH2D("hEtaPtRec", "hEtaPtRec", nEtaBin, etamin, etamax, nIncPtBin, logBinsX);
    hEtaERec = new TH2D("hEtaERec", "hEtaERec", nEtaBin, etamin, etamax, 200, 0., 200.);
    hEPhotonTrue = new TH1D("hEPhotonTrue", "hEPhotonTrue", 200, 0., 200.);
    hEPhotonCluster = new TH1D("hEPhotonCluster", "hEPhotonCluster", 200, 0., 200.);
}

void FillTruePions()
{
    int ntrack = fTracks->GetEntries();
    for (int itrack=0; itrack<ntrack; itrack++) {
        AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(itrack);
        double trEta = tr->Eta();
        double trPt = tr->Pt();
        double trE = tr->E();
        double trPid = tr->GetID();
        double trMomType = tr->GetMotherType();
        if (trPid==111) {
            hEtaPtTrue->Fill(trEta, trPt);
            hEtaETrue->Fill(trEta, trE);
        }
        if (trPid==22 && trMomType==111 && tr->IsPrimary())
            hEPhotonTrue->Fill(trE);
    }
}

void FillRecPions()
{
    int nclust = fTrackClusters->GetEntriesFast();
    for (int i = 1; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)fTrackClusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)fTrackClusters->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(fTrackClusters, lv1, lv2);
            if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())>asymcut) continue;
            double mass = 1000.*lvSum.M();
            if (mass > 110. && mass < 160.) {
                hEtaPtRec->Fill(lvSum.Eta(), lvSum.Pt());
                hEtaERec->Fill(lvSum.Eta(), lvSum.E());
            }
        }
    }
}

AliJBaseTrack GetPhotonSumVector(TClonesArray *clusters, AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    AliJBaseTrack lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}
