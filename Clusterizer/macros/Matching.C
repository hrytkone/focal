#include "Matching.h"

void Matching(TString inputfile)
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

            lvParticle.SetPtEtaPhiM(clustPt, clustEta, clustPhi, 0.);
            AliJBaseTrack track( lvParticle );
            new((*fTrackClusters)[fTrackClusters->GetEntriesFast()]) AliJBaseTrack(track);


            /**double dr = 10000.;
            int iClosestTrack = -1;

            int ntrack = fTracks->GetEntries();
            for (int itrack=0; itrack<ntrack; itrack++) {
                AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(itrack);
                double trEta = tr->Eta();
                double trPhi = tr->Phi();
                double trPt = tr->Pt();
                double trE = tr->E();
                double r = TMath::Sqrt((trPt-clustPt)*(trPt-clustPt) + (trE-clustE)*(trE-clustE) + (trEta-clustEta)*(trEta-clustEta) + (trPhi-clustPhi)*(trPhi-clustPhi));
                //if (r<dr && tr->GetMotherType()==111) {
                if (r<dr) {
                    dr = r;
                    iClosestTrack = itrack;
                }
            }**/
            /**if (iClosestTrack!=-1) {
                AliJBaseTrack *trClosest = (AliJBaseTrack*)fTracks->At(iClosestTrack);
                track_match.push_back(iClosestTrack);
                AliJBaseTrack *trClust = (AliJBaseTrack*)fTrackClusters->At(iclust);
                hPtECluster->Fill(clustPt, clustE);
                hEPerPtCluster->Fill(clustE/trClosest->Mag());

                double dz = fGeom->GetFOCALSegmentZ(1) - trClosest->Z();
                double x = trClosest->X() + trClosest->Px()/trClosest->Pz() * dz;
                double y = trClosest->Y() + trClosest->Py()/trClosest->Pz() * dz;
                //cout << "\nClosest track is " << iClosestTrack << " with PID=" << trClosest->GetID() << endl;
                //cout << "Cluster:\tpT=" << clustPt << "\teta=" << clustEta << "\tphi=" << clustPhi << "\t(x y)=(" << clustX << " " << clustY << ")\tE=" << clustE << endl;
                //cout << "Track:  \tpT=" << trClosest->Pt() << "\teta=" << trClosest->Eta() << "\tphi=" << trClosest->Phi() << "\t(x y)=(" << x << " " << y << ")\tE=" << trClosest->E() << "\tMotherPID=" << trClosest->GetMotherType() << "\tMohterID=" << trClosest->GetMotherID() << endl;
                if (trClosest->GetID()==11) {
                    trClust->SetID(11);
                    hEnergyElectron->Fill(trClosest->E(), clustE);
                    hPtElectron->Fill(trClosest->Pt(), clustPt);
                    hEnergyHCALElectron->Fill(clustEHCAL*0.000001);
                    hPtEElectron->Fill(clustPt, clustE);
                }
                if (TMath::Abs(trClosest->GetID())==211) {
                    trClust->SetID(211);
                    hEnergyChargedPion->Fill(trClosest->E(), clustE);
                    hPtChargedPion->Fill(trClosest->Pt(), clustPt);
                    hEnergyHCALChargedPion->Fill(clustEHCAL*0.000001);
                    hPtEChargedPion->Fill(clustPt, clustE);
                }
                if (trClosest->GetID()==22) {
                    trClust->SetID(22);
                    hEnergyPhoton->Fill(trClosest->E(), clustE);
                    hPtPhoton->Fill(trClosest->Pt(), clustPt);
                    hEnergyHCALPhoton->Fill(clustEHCAL*0.000001);
                    hPtEPhoton->Fill(clustPt, clustE);
                    if (trClosest->GetMotherType()==111) {
                        trClust->SetMotherType(111);
                        trClust->SetMotherID(trClosest->GetMotherID());
                        hEnergyPhotonFromDecay->Fill(trClosest->E(), clustE);
                        hPtPhotonFromDecay->Fill(trClosest->Pt(), clustPt);
                        hEnergyHCALPhotonFromDecay->Fill(clustEHCAL*0.000001);
                    }
                }
            } else {
                cout << "\nMATCH NOT FOUND" << endl;
                cout << "Cluster:\tpT=" << clustPt << "\teta=" << clustEta << "\tphi=" << clustPhi << "\tE=" << clustE << "\tHCAL=" << clustEHCAL << endl;
            }**/
        }
        FillMassHistos(fTrackClusters, fTracks);
        fTrackClusters->Clear("C");
        track_match.clear();
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
    hEnergyPhoton = new TH2D("hEnergyPhoton", "hEnergyPhoton", 200, 0., 100., 200, 0., 100.);
    hEnergyPhotonFromDecay = new TH2D("hEnergyPhotonFromDecay", "hEnergyPhotonFromDecay", 200, 0., 100., 200, 0., 100.);
    hEnergyChargedPion = new TH2D("hEnergyChargedPion", "hEnergyChargedPion", 200, 0., 100., 200, 0., 100.);
    hEnergyElectron = new TH2D("hEnergyElectron", "hEnergyElectron", 200, 0., 100., 200, 0., 100.);
    hPtPhoton = new TH2D("hPtPhoton", "hPtPhoton", 200, 0., 20., 200, 0., 20.);
    hPtPhotonFromDecay = new TH2D("hPtPhotonFromDecay", "hPtPhotonFromDecay", 200, 0., 20., 200, 0., 20.);
    hPtChargedPion = new TH2D("hPtChargedPion", "hPtChargedPion", 200, 0., 20., 200, 0., 20.);
    hPtElectron = new TH2D("hPtElectron", "hPtElectron", 200, 0., 20., 200, 0., 20.);

    hPtECluster = new TH2D("hPtECluster", "hPtECluster", 200, 0., 20., 200, 0., 100.);
    hPtEElectron = new TH2D("hPtEElectron", "hPtEElectron", 200, 0., 20., 200, 0., 100.);
    hPtEPhoton = new TH2D("hPtEPhoton", "hPtEPhoton", 200, 0., 20., 200, 0., 100.);
    hPtEChargedPion = new TH2D("hPtEChargedPion", "hPtEChargedPion", 200, 0., 20., 200, 0., 100.);

    hEPerPtCluster = new TH1D("hEPerPtCluster", "hEPerPtCluster", 100, 0., 5.);

    hEnergyHCALPhoton = new TH1D("hEnergyHCALPhoton", "hEnergyHCALPhoton", 250, 0., 50.);
    hEnergyHCALPhotonFromDecay = new TH1D("hEnergyHCALPhotonFromDecay", "hEnergyHCALPhotonFromDecay", 250, 0., 50.);
    hEnergyHCALChargedPion = new TH1D("hEnergyHCALChargedPion", "hEnergyHCALChargedPion", 250, 0., 50.);
    hEnergyHCALElectron = new TH1D("hEnergyHCALElectron", "hEnergyHCALElectron", 250, 0., 50.);

    for (int i=0; i<nasym; i++) {
        hCounter[i] = new TH1D(Form("hCounter_%d", i), "hCounter", npt, 0., npt);
        for (int j=0; j<npt; j++) {
            hMassCluster[i][j] = new TH1D(Form("hMassCluster_%d_%d", i, j), "hMassCluster", 350, 0., 700.);
            hMassPhoton[i][j] = new TH1D(Form("hMassPhoton_%d_%d", i, j), "hMassPhoton", 350, 0., 700.);
            hMassPhotonFromDecay[i][j] = new TH1D(Form("hMassPhotonFromDecay_%d_%d", i, j), "hMassPhotonFromDecay", 350, 0., 700.);
        }
    }
}

void FillMassHistos(TClonesArray *clusters, TClonesArray *tracks)
{
    int nclust = clusters->GetEntriesFast();
    for (int i = 1; i < nclust; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)clusters->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)clusters->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(clusters, lv1, lv2);
            if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())>asymcut[nasym-1]) continue; // to make
            if (lvSum.Eta()<4. || lvSum.Eta()>5.) continue;
            //if (lv1->E() < 5. || lv2->E() < 5.) continue;
            int ptbin = GetBin(pt, npt+1, lvSum.Pt());
            if (ptbin==-1) continue;
            for (int k=0; k<nasym; k++) {
                if (TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E())<asymcut[k]) {
                    double mass = 1000.*lvSum.M();
                    hMassCluster[k][ptbin]->Fill(mass);
                }
            }
            //if (lv1->GetID()==22 && lv2->GetID()==22) {
            //    hMassPhoton->Fill(mass);
            //    if (lv1->GetMotherType()==111 && lv2->GetMotherType()==111 && (lv1->GetMotherID()==lv2->GetMotherID())) {
            //        AliJBaseTrack *tr1 = (AliJBaseTrack*)tracks->At(track_match[i]);
            //        AliJBaseTrack *tr2 = (AliJBaseTrack*)tracks->At(track_match[j]);
            //        hMassPhotonFromDecay->Fill(mass);
            //        //cout << "\nMass=" << mass << "\tMass(calc)=" << 1000.*GetMass(lv1, lv2) << "\tMass(track)=" << 1000.*GetMass(tr1, tr2) << "\t" << tr1->GetMotherType() << " " << tr2->GetMotherType() << "\tpT=" << lvSum.Pt() <<  endl;
            //        //cout << "Track 1:  \tpT=" << tr1->Pt() << "\teta=" << tr1->Eta() << "\tphi=" << tr1->Phi() << "\tE=" << tr1->E() << endl;
            //        //cout << "Cluster 1:  \tpT=" << lv1->Pt() << "\teta=" << lv1->Eta() << "\tphi=" << lv1->Phi() << "\tE=" << lv1->E() << endl;
            //        //cout << "Track 2:  \tpT=" << tr2->Pt() << "\teta=" << tr2->Eta() << "\tphi=" << tr2->Phi() << "\tE=" << tr2->E() << endl;
            //        //cout << "Cluster 2:  \tpT=" << lv2->Pt() << "\teta=" << lv2->Eta() << "\tphi=" << lv2->Phi() << "\tE=" << lv2->E() << endl;
            //    }
            //}
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

AliJBaseTrack GetPhotonSumVector(TClonesArray *clusters, AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    AliJBaseTrack lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}

double GetMass(AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    TVector3 v1(lv1->Px(), lv1->Py(), lv1->Pz());
    TVector3 v2(lv2->Px(), lv2->Py(), lv2->Pz());
    double angle = TMath::ACos((v1.Dot(v2))/(v1.Mag()*v2.Mag()));
    return TMath::Sqrt(2.*lv1->E()*lv2->E()*(1.-TMath::Cos(angle)));
}
