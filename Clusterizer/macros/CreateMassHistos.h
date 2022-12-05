//const TString outputname = "test";
//const TString outputname = "etacut_37-56_pthard-2";
//const TString outputname = "etacut_41-55_more-pt-bins";
const TString outputname = "etacut_41-55_leading-trigg-eff_1";
const int nasym = 6;
//const int npt = 6;
const int npt = 10;

const double etamin = 4.1;
const double etamax = 5.5;
//const double etamin = 4.;
//const double etamax = 5.;

const int poolsize = 0;

double asymcut[nasym] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
//double pt[npt+1] = {1., 2., 3., 5., 8., 12., 20.};
//double pt[npt+1] = {2., 3., 4., 8., 10., 15., 20.};
double pt[npt+1] = {2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 15., 20.};
//double pt[npt+1] = {10., 11.};

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;

AliFOCALGeometry *fGeom;

TClonesArray *fTracks;
TClonesArray *fClusters;
TClonesArray *fTrackClusters;
TClonesArray *fClustersMixed[poolsize];
TLorentzVector lvParticle;

// Histograms
TH1D *hCounter;
TH1D *hCounterTrue;
TH1D *hIncPt;
TH1D *hMassCluster[nasym][npt];
TH1D *hMassClusterMixed[nasym][npt];
TH1D *hMassClusterRotated[nasym][npt];

Int_t LoadInput(TString inputfile);
void InitOutput();
void FillMassHistos(TClonesArray *clusters, TClonesArray *tracks);
void FillMassHistosMixed(TClonesArray *clusters, TClonesArray *clustpool);
void FillMassHistosRotated(TClonesArray *clusters);
AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2);
int GetBin(double arr[], int nArr, double val);