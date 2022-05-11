const TString outputname = "eta-cut";
const int nasym = 6;
const int npt = 6;

//double ecut = 2.;

double asymcut[nasym] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
double pt[npt+1] = {1., 2., 3., 5., 8., 12., 20.};

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;

AliFOCALGeometry *fGeom;

TClonesArray *fTracks;
TClonesArray *fClusters;
TClonesArray *fTrackClusters = new TClonesArray("AliJBaseTrack", 1500);
TLorentzVector lvParticle;

// Histograms
TH2D *hEnergyPhoton;
TH2D *hEnergyPhotonFromDecay;
TH2D *hEnergyElectron;
TH2D *hEnergyChargedPion;
TH2D *hPtPhoton;
TH2D *hPtPhotonFromDecay;
TH2D *hPtElectron;
TH2D *hPtChargedPion;
TH2D *hPtECluster;
TH2D *hPtEElectron;
TH2D *hPtEPhoton;
TH2D *hPtEChargedPion;

TH1D *hEPerPtCluster;

TH1D *hEnergyHCALPhoton;
TH1D *hEnergyHCALPhotonFromDecay;
TH1D *hEnergyHCALElectron;
TH1D *hEnergyHCALChargedPion;

TH1D *hCounter[nasym][npt];
TH1D *hMassCluster[nasym][npt];
TH1D *hMassPhoton[nasym][npt];
TH1D *hMassPhotonFromDecay[nasym][npt];

vector<int> track_match;

Int_t LoadInput(TString inputfile);
void InitOutput();
void FillMassHistos(TClonesArray *clusters, TClonesArray *tracks);
AliJBaseTrack GetPhotonSumVector(TClonesArray *clusters, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
double GetMass(AliJBaseTrack *lv1, AliJBaseTrack *lv2);
int GetBin(double arr[], int nArr, double val);
