void FillMassClusterFiles(TString inputdir, TString tag)
{
    for (int ifile=0; ifile<nfile; ifile++) {
        TFile *fIn = TFile::Open(Form("%s/%s_%d-%d/clusters_MBPythia_%d.root", inputdir.Data(), tag.Data(), ifile, ifile, ifile), "READ");
        if (!fIn) {
            cout << "File not found!" << endl;
            break;
        }

        TTree *fTree = (TTree*)fIn->Get("Run");
        fTree->SetBranchAddress("Events", &fEvent);
        fBrtrack = fTree->GetBranch("Events.fTracks");
        fBrtrack->SetAddress(&fTracks);
        fBrcluster = fTree->GetBranch("Events.fClusters");
        fBrcluster->SetAddress(&fClusters);
    }
}
//******************************************************************************
