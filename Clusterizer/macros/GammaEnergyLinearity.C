TString outputname = "gamma-energy_12900";

TTree* fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
TFile *fIn, *fOut;
TClonesArray *fTracks;
TClonesArray *fClusters;
AliJHMREvent *fEvent;

TH2D *hEnergy;

Int_t LoadInput(TString inputfile);
void InitOutput();

void GammaEnergyLinearity(TString inputfile)
{
    if (!LoadInput(inputfile)) {
        cout << "File not found! " << inputfile << endl;
        return;
    }

    InitOutput();

    int nev = fTree->GetEntries();
    cout << "Processing " << nev << " events" << endl;
    for (int iev=0; iev<nev; iev++) {

        if (iev%100==0) cout << "event " << iev << "/" << nev << endl;
        //cout << "\n================================" << endl;

        fTree->GetEntry(iev);

        float ecluster = 0;

        int nclust = fClusters->GetEntries();
        for (int iclust=0; iclust<nclust; iclust++) {
            AliJHMRCluster *clust = (AliJHMRCluster*)fClusters->At(iclust);
            ecluster += clust->GetE();
        }

        AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(0);
        double emc = tr->E();

        hEnergy->Fill(ecluster, emc);
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
    hEnergy = new TH2D("hEnergy", "hEnergy", 160, 0., 1600., 160, 0., 1600.); hEnergy->Sumw2();
}
