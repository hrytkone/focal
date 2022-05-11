//const TString outputname = "numberOfClusters_lclust-ecut-3";
const TString outputname = "numberOfClusters";
const int npt = 6;
const int necut = 6;
const int nptbin = 120;

double pt[npt+1] = {2., 3., 4., 8., 10., 15., 20.};
double ecut[necut] = {0., 0.1, 0.35, 0.5, 1., 3.};

double logBinsX[nptbin+1], limMin = 0.1, limMax = 20;
const double logBW = (log(limMax) - log(limMin))/nptbin;

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;

TClonesArray *fTracks;
TClonesArray *fClusters;

// Histograms
TH1D *hNumberOfClusters[npt][necut];
TH2D *hClusterPerGammaEnergy;
TH2D *hLeadingClusterDeltaPhiDeltaEta;
TH1D *hClusterMoreThanOne;
TH1D *hClusterPt;

Int_t LoadInput(TString inputfile);
void InitOutput();
int GetBin(double arr[], int nArr, double val);

void CheckNumberOfClustersForGammas(TString inputfile)
{

    if (!LoadInput(inputfile)) {
        cout << "File not found! " << inputfile << endl;
        return;
    }

    InitOutput();

    int nev = fTree->GetEntries();
    cout << "Processing " << nev << " events" << endl;
    for (int iev=0; iev<nev; iev++) {

        //if (iev%100==0) cout << "event " << iev << "/" << nev << endl;
        //cout << "\n================================" << endl;

        fTree->GetEntry(iev);

        int ntrack = fTracks->GetEntries();
        AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(0);

        double leadingClustE = -100.;
        double leadingClustPhi = -100.;
        double leadingClustEta = -100.;
        int nclustAfterCut[necut] = {0};
        int nclust = fClusters->GetEntries();
        //cout << "\nGamma energy = " << tr->E() << "\tnclust = " << nclust << endl;
        for (int iclust=0; iclust<nclust; iclust++) {
            AliJHMRCluster *clust = (AliJHMRCluster*)fClusters->At(iclust);
            float clustE = clust->GetE();
            float clustEta = clust->GetEta();
            float clustPhi = clust->GetPhi();
            //cout << "\tE: " << clustE << endl;
            if (clustE > leadingClustE) {
                leadingClustE = clustE;
                leadingClustEta = clustEta;
                leadingClustPhi = clustPhi;
            }
            for (int iecut=0; iecut<necut; iecut++) {
                if (clustE > ecut[iecut]) nclustAfterCut[iecut]++;
            }
        }

        double trPt = tr->Pt();
        double trE = tr->E();
        double trEta = tr->Eta();
        double trPhi = tr->Phi();
        int ibin = GetBin(pt, npt, trPt);
        if (ibin > -1) {
            for (int iecut=0; iecut<necut; iecut++) {
                hNumberOfClusters[ibin][iecut]->Fill(nclustAfterCut[iecut]);
            }
            hClusterPt->Fill(trPt);
            if (nclust>1) hClusterMoreThanOne->Fill(trPt);
        }

        // Save ratio of leading cluster energy to gamma energy
        double ratio = leadingClustE/trE;
        //if (leadingClustE<3.0) continue;
        hClusterPerGammaEnergy->Fill(trPt, ratio);
        hLeadingClusterDeltaPhiDeltaEta->Fill(leadingClustPhi-trPhi, leadingClustEta-trEta);
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
    for (int i=0; i<npt; i++) {
        for (int j=0; j<necut; j++) {
            hNumberOfClusters[i][j] = new TH1D(Form("hNumberOfClusters_%d_%d", i,j), "hNumberOfClusters", 30, -0.5, 29.5); hNumberOfClusters[i][j]->Sumw2();
        }
    }
    for (int i = 0; i <= nptbin; i++) logBinsX[i] = limMin*exp(i*logBW);
    hClusterPerGammaEnergy = new TH2D("hClusterPerGammaEnergy", "", nptbin, logBinsX, nptbin, 0., 1.2); hClusterPerGammaEnergy->Sumw2();
    hLeadingClusterDeltaPhiDeltaEta = new TH2D("hLeadingClusterDeltaPhiDeltaEta", "", 250, -TMath::Pi(), TMath::Pi(), 250, -2.5, 2.5); hLeadingClusterDeltaPhiDeltaEta->Sumw2();
    hClusterPt = new TH1D("hClusterPt", "hClusterPt", npt, pt); hClusterPt->Sumw2();
    hClusterMoreThanOne = new TH1D("hClusterMoreThanOne", "hClusterMoreThanOne", npt, pt); hClusterMoreThanOne->Sumw2();
}

int GetBin(double arr[], int nArr, double val)
{
    for (int i=0; i<nArr; i++) {
        if (arr[i]<=val && val<arr[i+1]) return i;
    }
    return -1;
}
