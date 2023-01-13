#include "TriggerPtModification.h"

void TriggerPtModification(TString inputfile)
{
    fGeom = AliFOCALGeometry::GetInstance("../geometry.txt");

    if (!LoadInput(inputfile)) {
        cout << "File not found! " << inputfile << endl;
        return;
    }

    InitOutput();

    int nev = fTree->GetEntries();
    //nev = 100000;
    cout << "Processing " << nev << " events" << endl;
    for (int iev=0; iev<nev; iev++) {
        if (iev%1000==0) cout << "event " << iev << "/" << nev << endl;

        fTree->GetEntry(iev);

        GetTruePions();
        FillTruePt();
        GetRecPions();
        FillRecPt();

        fTruePi0->Clear("C");
        fRecPi0->Clear("C");
    }
    fOut->cd();
    fOut->Write("", TObject::kOverwrite);
}

//************************************************************************************************
//************************************************************************************************

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
    hPtTrue = new TH1D("hPtTrue", "hPtTrue", nPtBin, ptMin, ptMax); hPtTrue->Sumw2();
    hPtRec = new TH1D("hPtRec", "hPtRec", nPtBin, ptMin, ptMax); hPtRec->Sumw2();
}

void GetTruePions()
{
    for (int itrack = 0; itrack < fTracks->GetEntries(); itrack++) {
        AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(itrack);
        double trEta = tr->Eta();
        // Use smaller acceptance than for gammas to suppress the effect from missing gamma pairs
        if ( trEta < etaMin+etacut || trEta > etaMax-etacut ) continue;

        if ( tr->GetID() == 111 ) {
            AliJBaseTrack track = *tr;
            new((*fTruePi0)[fTruePi0->GetEntriesFast()]) AliJBaseTrack(track);
        }
    }
}

void GetRecPions()
{
    TLorentzVector lv1, lv2;
    for (int i = 1; i < fClusters->GetEntriesFast(); i++) {
        AliJHMRCluster *c1 = (AliJHMRCluster*)fClusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJHMRCluster *c2 = (AliJHMRCluster*)fClusters->At(j);
            if (GetAsymmetry(c1->GetE(), c2->GetE())>asymcut) continue;
            lv1.SetPtEtaPhiM(c1->GetPt(), c1->GetEta(), c1->GetPhi(), 0.);
            lv2.SetPtEtaPhiM(c2->GetPt(), c2->GetEta(), c2->GetPhi(), 0.);
            AliJBaseTrack tr1(lv1);
            AliJBaseTrack tr2(lv2);
            AliJBaseTrack lvSum = GetPhotonSumVector(&tr1, &tr2);
            if (lvSum.Eta()<etaMin+etacut || lvSum.Eta()>etaMax-etacut) continue;
            double mass = 1000.*lvSum.M();
            if (mass < mwMin || mass > mwMax) continue;
            new ((*fRecPi0)[fRecPi0->GetEntriesFast()]) AliJBaseTrack(lvSum);
        }
    }
}


void FillTruePt()
{
    for (int itrack=0; itrack<fTruePi0->GetEntries(); itrack++) {
        AliJBaseTrack *tr = (AliJBaseTrack*)fTruePi0->At(itrack);
        double trPt = tr->Pt();
        hPtTrue->Fill(trPt);
    }
}

void FillRecPt()
{
    for (int iclust=0; iclust<fRecPi0->GetEntries(); iclust++) {
        AliJBaseTrack *clust = (AliJBaseTrack*)fRecPi0->At(iclust);
        double clustPt = clust->Pt();
        hPtRec->Fill(clustPt);
    }
}

AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    AliJBaseTrack lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}

double GetAsymmetry(double E1, double E2)
{
    return TMath::Abs(E1 - E2)/(E1 + E2);
}