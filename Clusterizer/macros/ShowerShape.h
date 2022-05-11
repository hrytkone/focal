const TString outputname = "shower-shape_gamma";
const int nseg = 10;

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;

TClonesArray *fTracks;
TClonesArray *fClusters;

TH1D *hWidth1[nseg];
TH1D *hWidth2[nseg];
TH1D *hWidth2PerWidth1[nseg];

Int_t LoadInput(TString inputfile);
void InitOutput();
