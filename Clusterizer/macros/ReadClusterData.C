TFile *fOut;
TH1D* hMassCluster;

double asymcut = 0.8;

void InitOutput();
TLorentzVector GetPhotonSumVector(TLorentzVector *lv1, TLorentzVector *lv2);
void FillRecPions(TClonesArray *clusters);

void ReadClusterData(TString dir, TString tag, int startFolder, int endFolder, int startEvent, int endEvent)
{
    AliFOCALGeometry * geometry = AliFOCALGeometry::GetInstance("../geometry.txt");

    InitOutput();

    // Loop over folders
    for (Int_t ifolder = startFolder; ifolder <= endFolder; ifolder++) {

        TString filename = Form("%s/%s_%d-%d/clusters_MBPythia_%d.root", dir.Data(), tag.Data(), ifolder, ifolder, ifolder);
        cout << "FOLDER: " << ifolder << "\tFILE: " << filename << endl;

        TFile *inputFile = new TFile(filename.Data(), "READ");
        TTree * tClusters = 0;

        // LOOP OVER EVENTS
        for (Int_t ievt = startEvent; ievt <= endEvent; ievt++) {

            inputFile->GetDirectory(Form("Event%i",ievt))->GetObject("fTreeR", tClusters);

            TBranch * bClusters = tClusters->GetBranch("AliFOCALCluster");

            TClonesArray * clustersArray = 0;
            bClusters->SetAddress(&clustersArray);
            bClusters->GetEvent(0);

        } // END OF EVENT LOOP

    }
}

void InitOutput()
{
    fOut = TFile::Open("output_clusters.root", "RECREATE");
    hMassCluster = new TH1D("hMassCluster", "hMassCluster", 400, 0., 400.);
}

TLorentzVector GetPhotonSumVector(TLorentzVector *lv1, TLorentzVector *lv2)
{
    TLorentzVector lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}

void FillRecPions(TClonesArray *clusters)
{
    int nclust = clusters->GetEntriesFast();
    for (int i=1; i<nclust; i++) {
        AliFOCALCluster *cluster1 = (AliFOCALCluster*) clusters->UncheckedAt(i);
        Double_t e1 = cluster1->E();
        Double_t x1 = cluster1->X();
        Double_t y1 = cluster1->Y();
        Double_t z1 = cluster1->Z();
        Double_t d1 = TMath::Sqrt(x1*x1 + y1*y1 + z1*z1);
        TLorentzVector lv1;
        lv1.SetPxPyPzE(e1*x1/d1, e1*y1/d1, e1*z1/d1, e1);
        for (int j=0; j<i; j++) {
            AliFOCALCluster *cluster2 = (AliFOCALCluster*) clusters->UncheckedAt(i);
            Double_t e2 = cluster2->E();
            Double_t x2 = cluster2->X();
            Double_t y2 = cluster2->Y();
            Double_t z2 = cluster2->Z();
            Double_t d2 = TMath::Sqrt(x2*x2 + y2*y2 + z2*z2);
            TLorentzVector lv2;
            lv2.SetPxPyPzE(e2*x2/d2, e2*y2/d2, e2*z2/d2, e2);

            AliJBaseTrack lvSum = GetPhotonSumVector(lv1, lv2);
            if (TMath::Abs(e1 - e2)/(e1 + e2) < asymcut) {
                double mass = 1000.*lvSum.M();
                hMassCluster->Fill(mass);
            }
        }
    }
}
