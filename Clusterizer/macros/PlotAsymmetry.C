const TString outputname = "asymmetry.root";

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;
TClonesArray *fTracks;
TClonesArray *fClusters;

TH1D *hAsymTrue;
TH1D *hAsymClust;

Int_t LoadInput(TString inputfile);
void InitOutput();

void PlotAsymmetry(TString inputfile)
{
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
        for (int i=1; i<nclust; i++) {
            AliJHMRCluster *clust1 = (AliJHMRCluster*)fClusters->At(i);
            for (int j=0; j<i; j++) {
                AliJHMRCluster *clust2 = (AliJHMRCluster*)fClusters->At(j);
                double e1 = clust1->GetE();
                double e2 = clust2->GetE();
                hAsymClust->Fill(TMath::Abs(e1 - e2)/(e1 + e2));
            }
        }

        int ntrack = fTracks->GetEntries();
        for (int i=1; i<ntrack; i++) {
            AliJBaseTrack *tr1 = (AliJBaseTrack*)fTracks->At(i);
            int imom1 = tr1->GetMotherID();
            if (imom1!=111) continue;
            if (!tr1->IsPrimary()) continue;
            for (int j=0; j<i; j++) {
                AliJBaseTrack *tr2 = (AliJBaseTrack*)fTracks->At(j);
                double e1 = tr1->E();
                double e2 = tr2->E();
                int imom2 = tr2->GetMotherID();
                if (imom2!=111) continue;
                if (!tr2->IsPrimary()) continue;
                hAsymTrue->Fill(TMath::Abs(e1 - e2)/(e1 + e2));
            }
        }
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
    hAsymClust = new TH1D("hAsymClust", "hAsymClust", 100, 0., 1.);
    hAsymTrue = new TH1D("hAsymTrue", "hAsymTrue", 100, 0., 1.);
}
