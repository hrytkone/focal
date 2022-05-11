#include "CreateMassHistos.h"

void CreateMassHistos(TString inputfile)
{
    fGeom = AliFOCALGeometry::GetInstance("../geometry.txt");

    if (!LoadInput(inputfile)) {
        cout << "File not found! " << inputfile << endl;
        return;
    }

    InitOutput();

    //int nev = 1;
    int nev = fTree->GetEntries();
    cout << "Processing " << nev << " events" << endl;
    for (int iev=0; iev<nev; iev++) {

        if (iev%100==0) cout << "event " << iev << "/" << nev << endl;
        //cout << "\n================================" << endl;

        fTree->GetEntry(iev);

        int nclust = fClusters->GetEntries();
        for (int iclust=0; iclust<nclust; iclust++) {
            AliJHMRCluster *clust = (AliJHMRCluster*)fClusters->At(iclust);
            float clustEta = clust->GetEta();
            float clustPhi = clust->GetPhi();
            float clustX = clust->GetX();
            float clustY = clust->GetY();
            float clustE = clust->GetE();
            float clustEHCAL = clust->GetEHCAL();
            float clustPt = clust->GetPt();
            float w1 = clust->GetWidth1(0);
            float w2 = clust->GetWidth2(0);

            hIncPt->Fill(clustPt);

            lvParticle.SetPtEtaPhiM(clustPt, clustEta, clustPhi, 0.);
            AliJBaseTrack track( lvParticle );
            new((*fTrackClusters)[fTrackClusters->GetEntriesFast()]) AliJBaseTrack(track);
        }

        //vector<int> motherID;
        int ntrack = fTracks->GetEntries();
        for (int itrack=0; itrack<ntrack; itrack++) {
            AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(itrack);
            double trEta = tr->Eta();
            double trPt = tr->Pt();
            double trPid = tr->GetID();
            if (trPid==111 && (trEta > etamin && trEta < etamax)) {
                int ptbin = GetBin(pt, npt+1, trPt);
                hCounter->Fill(ptbin+0.5);
            }
            /**double trMotherID = tr->GetMotherID();
            if (trMotherID<0) continue;
            if (trPid==22 && tr->IsPrimary() && (trEta > etamin && trEta < etamax)) {
                if (std::count(motherID.begin(), motherID.end(), trMotherID)) {
                    AliJBaseTrack *trMother = (AliJBaseTrack*)fTracks->At(trMotherID);
                    int ptbin = GetBin(pt, npt+1, trMother->Pt());
                    hCounterTrue->Fill(ptbin+0.5);
                } else {
                    motherID.push_back(trMotherID);
                }
            }**/
        }
        //motherID.clear();

        FillMassHistos(fTrackClusters, fTracks);
        //FillMassHistosRotated(fTrackClusters);

        // Mixed event : take triggers from this event, associated from previous
        if (poolsize>0) {
            if (iev >= poolsize) {
                for (int ipool = 0; ipool < poolsize; ipool++) {
                    FillMassHistosMixed(fTrackClusters, fClustersMixed[ipool]);
                }
                // Remove the first from the pool and add new array to the pool
                for (int ipool = 0; ipool < poolsize-1; ipool++) {
                    *fClustersMixed[ipool] = *fClustersMixed[ipool+1];
                }
                *fClustersMixed[poolsize-1] = *fTrackClusters;
            } else { // Create pools
                *fClustersMixed[iev] = *fTrackClusters;
            }
        }

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
    hCounter = new TH1D("hCounter", "hCounter", npt, 0., npt); hCounter->Sumw2();
    hCounterTrue = new TH1D("hCounterTrue", "hCounterTrue", npt, 0., npt); hCounterTrue->Sumw2();
    hIncPt = new TH1D("hIncPt", "hIncPt", 100, 0., 20.); hIncPt->Sumw2();
    for (int i=0; i<nasym; i++) {
        for (int j=0; j<npt; j++) {
            hMassCluster[i][j] = new TH1D(Form("hMassCluster_%d_%d", i, j), "hMassCluster", 350, 0., 700.); hMassCluster[i][j]->Sumw2();
            hMassClusterMixed[i][j] = new TH1D(Form("hMassClusterMixed_%d_%d", i, j), "hMassClusterMixed", 350, 0., 700.); hMassClusterMixed[i][j]->Sumw2();
            hMassClusterRotated[i][j] = new TH1D(Form("hMassClusterRotated_%d_%d", i, j), "hMassClusterRotated", 350, 0., 700.); hMassClusterRotated[i][j]->Sumw2();
        }
    }

    fTrackClusters = new TClonesArray("AliJBaseTrack", 1500);
    for (int ipool = 0; ipool < poolsize; ipool++)
        fClustersMixed[ipool] = new TClonesArray("AliJBaseTrack", 1500);
}

void FillMassHistos(TClonesArray *clusters, TClonesArray *tracks)
{
    int nclust = clusters->GetEntriesFast();
    for (int i = 1; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)clusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)clusters->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(lv1, lv2);
        //    if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())>asymcut[nasym-1]) continue;
            if (lvSum.Eta()<etamin || lvSum.Eta()>etamax) continue;
            int ptbin = GetBin(pt, npt+1, lvSum.Pt());
            if (ptbin==-1) continue;
            for (int k=0; k<nasym; k++) {
                if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())<asymcut[k]) {
                    double mass = 1000.*lvSum.M();
                    hMassCluster[k][ptbin]->Fill(mass);
                }
            }
        }
    }
}

void FillMassHistosRotated(TClonesArray *clusters)
{
    int nclust = clusters->GetEntriesFast();
    for (int i = 1; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)clusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)clusters->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(lv1, lv2);

            AliJBaseTrack lv1rot = *lv1;
            AliJBaseTrack lv2rot = *lv2;

            lv1rot.Rotate(TMath::Pi()/2., lvSum.Vect());
            lv2rot.Rotate(TMath::Pi()/2., lvSum.Vect());

            for (int k = 0; k < nclust; k++) {
                if (k==i || k==j) continue;
                AliJBaseTrack *lv3 = (AliJBaseTrack*)clusters->At(k);
                AliJBaseTrack lvSumRot1 = GetPhotonSumVector(&lv1rot, lv3);
                if (lvSumRot1.Eta()>etamin && lvSumRot1.Eta()<etamax) {
                    int ptbin = GetBin(pt, npt+1, lvSumRot1.Pt());
                    if (ptbin!=-1) {
                        for (int l=0; l<nasym; l++) {
                            if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())<asymcut[l]) {
                                double mass = 1000.*lvSumRot1.M();
                                hMassClusterRotated[l][ptbin]->Fill(mass);
                            }
                        }
                    }

                }
                AliJBaseTrack lvSumRot2 = GetPhotonSumVector(&lv2rot, lv3);
                if (lvSumRot2.Eta()>etamin && lvSumRot2.Eta()<etamax) {
                    int ptbin = GetBin(pt, npt+1, lvSumRot2.Pt());
                    if (ptbin!=-1) {
                        for (int l=0; l<nasym; l++) {
                            if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())<asymcut[l]) {
                                double mass = 1000.*lvSumRot2.M();
                                hMassClusterRotated[l][ptbin]->Fill(mass);
                            }
                        }
                    }
                }
            }
        }
    }
}

void FillMassHistosMixed(TClonesArray *clusters, TClonesArray *clustpool)
{
    int nclust = clusters->GetEntriesFast();
    int nclustpool = clustpool->GetEntriesFast();
    for (int i = 0; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)clusters->At(i);
        for (int j = 0; j < nclustpool; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)clustpool->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(lv1, lv2);
            if (lvSum.Eta()<etamin || lvSum.Eta()>etamax) continue;
            int ptbin = GetBin(pt, npt+1, lvSum.Pt());
            if (ptbin==-1) continue;
            for (int k=0; k<nasym; k++) {
                if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())<asymcut[k]) {
                    double mass = 1000.*lvSum.M();
                    hMassClusterMixed[k][ptbin]->Fill(mass);
                }
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

AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    AliJBaseTrack lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}
