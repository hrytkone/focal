const TString outputname = "masses_v9";
const int nasym = 6;
const int npt = 4;

const int poolsize = 0;

const double etamin = 3.5;
const double etamax = 5.5;

const double etacut = 0.2;

double asymcut[nasym] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
double pt[npt+1] = {2., 3., 4., 8., 20.};

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
