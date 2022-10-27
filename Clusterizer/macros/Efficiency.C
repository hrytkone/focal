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
        FillMassHistos(fTrackClusters);
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

    for (int i=0; i<nEtaBin; i++) {
        for (int j=0; j<nPtBin; j++) {
            hMassCluster[i][j] = new TH1D(Form("hMassCluster_%d_%d", i, j), "hMassCluster", 350, 0., 700.); hMassCluster[i][j]->Sumw2();
        }
    }

    hEtaPtTrue = new TH2D("hEtaPtTrue", "hEtaPtTrue", nEtaBin, eta, nPtBin, pt);
    hEtaETrue = new TH2D("hEtaETrue", "hEtaETrue", nEtaBin, eta, 200, 0., 200.);
    hEtaPtRec = new TH2D("hEtaPtRec", "hEtaPtRec", nEtaBin, eta, nPtBin, pt);
    hEtaERec = new TH2D("hEtaERec", "hEtaERec", nEtaBin, eta, 200, 0., 200.);
    hEPhotonTrue = new TH1D("hEPhotonTrue", "hEPhotonTrue", 200, 0., 200.);
    hEPhotonCluster = new TH1D("hEPhotonCluster", "hEPhotonCluster", 200, 0., 200.);

    hPtMass  = new TH2D("hPtMass", "hPtMass", 36, 2., 20., 350, 0., 700.);
    for (int ipt=0; ipt<nPtBin; ipt++)
        hEtaMass[ipt] = new TH2D(Form("hEtaMass_%d", ipt), "hEtaMass", nEtaBin, eta, 350, 0., 700.);

    hPhiEtaTrue = new TH2D("hPhiEtaTrue", "hPhiEtaTrue", nPhiBin, phimin, phimax, nEtaBin, eta);
    hPhiEta = new TH2D("hPhiEta", "hPhiEta", nPhiBin, phimin, phimax, nEtaBin, eta);
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
        if (trPid==111) {
            hEtaPtTrue->Fill(trEta, trPt);
            hEtaETrue->Fill(trEta, trE);
            hPhiEtaTrue->Fill(trPhi, trEta);
        }
        if (trPid==22 && trMomType==111 && tr->IsPrimary())
            hEPhotonTrue->Fill(trE);
    }
}

void FillMassHistos(TClonesArray *clusters)
{
    int nclust = clusters->GetEntriesFast();
    for (int i = 1; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)clusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)clusters->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(lv1, lv2);
        //    if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())>asymcut[nasym-1]) continue;
            if (lvSum.Eta()<etamin || lvSum.Eta()>etamax) continue;
            int ptbin = GetBin(pt, nPtBin+1, lvSum.Pt());
            int etabin = GetBin(eta, nEtaBin+1, lvSum.Eta());
            if (ptbin==-1 || etabin==-1) continue;
            if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())<asymcut) {
                double mass = 1000.*lvSum.M();
                hMassCluster[etabin][ptbin]->Fill(mass);
                hPtMass->Fill(lvSum.Pt(), mass);
                hEtaMass[ptbin]->Fill(lvSum.Eta(), mass);
                if (mass > 110. && mass < 140.)
                    hPhiEta->Fill(lvSum.Phi(), lvSum.Eta());
            }
        }
    }
}

int GetBin(double arr[], int nArr, double val)
{
    for (int i=0; i<nArr; i++) {
        if (arr[i]<=val && val<arr[i+1]) return i;
    }
    return -1;
}

void FillRecPions()
{
    int nclust = fTrackClusters->GetEntriesFast();
    for (int i = 1; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)fTrackClusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)fTrackClusters->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(lv1, lv2);
            if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())>asymcut) continue;
            double mass = 1000.*lvSum.M();
            if (mass > 110. && mass < 160.) {
                hEtaPtRec->Fill(lvSum.Eta(), lvSum.Pt());
                hEtaERec->Fill(lvSum.Eta(), lvSum.E());
            }
        }
    }
}

AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    AliJBaseTrack lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}
