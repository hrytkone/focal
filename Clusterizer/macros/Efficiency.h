const TString outputname = "efficiency_asym-08_v1.3_calib";
const double asymcut = 0.8;

const int nPtBin = 6;
double pt[nPtBin+1], limMin = 2, limMax = 20;
double logBW = (log(limMax) - log(limMin))/nPtBin;

const int nEtaBin = 52;
double eta[nEtaBin+1];
double etaBW = 0.05, etamin = 3.2, etamax = 5.8;

const int nPhiBin = 52;
double phimin = -TMath::Pi(), phimax = TMath::Pi();

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
TH2D *hEtaPtTrue;
TH2D *hEtaPtRec;
TH2D *hEtaETrue;
TH2D *hEtaERec;
TH2D *hPtMass;
TH2D *hEtaMass[nPtBin];
TH2D *hPhiEtaTrue;
TH2D *hPhiEta;

TH1D *hMassCluster[nEtaBin][nPtBin];
TH1D *hEPhotonCluster;
TH1D *hEPhotonTrue;

Int_t LoadInput(TString inputfile);
void InitOutput();
void FillTruePions();
void FillRecPions();
void FillMassHistos(TClonesArray *clusters);
int GetBin(double arr[], int nArr, double val);
AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2);
