#include "ShowerShape.h"

void ShowerShape(TString inputfile)
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
        for (int iclust=0; iclust<nclust; iclust++) {
            AliJHMRCluster *clust = (AliJHMRCluster*)fClusters->At(iclust);
            for (int iseg = 0; iseg < nseg; iseg++) {
                float width1 = clust->GetWidth1(iseg);
                float width2 = clust->GetWidth2(iseg);
                if (width1 <= 0 || width2 <= 0) continue;
                hWidth1[iseg]->Fill(width1);
                hWidth2[iseg]->Fill(width2);
                hWidth2PerWidth1[iseg]->Fill(width2/width1);
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
    for (int iseg=0; iseg<nseg; iseg++) {
        hWidth1[iseg] = new TH1D(Form("hWidth1_%d", iseg), "", 250, 0., 5.);
        hWidth2[iseg] = new TH1D(Form("hWidth2_%d", iseg), "", 250, 0., 5.);
        hWidth2PerWidth1[iseg] = new TH1D(Form("hWidth2PerWidth1_%d", iseg), "", 200, 0., 1.);
    }
}
