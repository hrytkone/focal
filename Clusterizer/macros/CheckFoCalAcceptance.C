#include "CheckFoCalAcceptance.h"

void CheckFoCalAcceptance(TString inputfile)
{
    fGeom = AliFOCALGeometry::GetInstance("../geometry.txt");

    if (!LoadInput(inputfile)) {
        cout << "File not found! " << inputfile << endl;
        return;
    }

    InitOutput();

    int nev = 10000;
    //int nev = fTree->GetEntries();
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

            if (clustY>ycutmin && clustY<ycutmax)
                hXEta->Fill(clustX, clustEta);
            if (clustX>xcutmin && clustX<xcutmax)
                hYEta->Fill(clustY, clustEta);
            hXY->Fill(clustX, clustY);
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
    for (int i = 0; i <= nIncPtBin; i++) logBinsX[i] = limMin*exp(i*logBW);
    hXEta = new TH2D("hXEta", "hXEta", nXBin, xmin, xmax, nEtaBin, etamin, etamax);
    hYEta = new TH2D("hYEta", "hYEta", nXBin, xmin, xmax, nEtaBin, etamin, etamax);
    hXY = new TH2D("hXY", "hXY", nXBin, xmin, xmax, nXBin, xmin, xmax);
}
