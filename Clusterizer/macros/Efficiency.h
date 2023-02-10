const TString outputname = "efficiency_pi0-gun_MW-100-150_matched";
const double asymcut = 1.;

const int nPtBin = 38;
double pt[nPtBin+1], limMin = 2, limMax = 18;
double logBW = (log(limMax) - log(limMin))/nPtBin;

const int nEtaBin = 38;
double eta[nEtaBin+1], etamin = 3.4, etamax = 5.3;
double etaBW = (etamax - etamin)/nEtaBin;

//const int nPhiBin = 52;
const int nPhiBin = 104;
double phimin = -TMath::Pi(), phimax = TMath::Pi();

//const int nThetaBin = 52;
const int nThetaBin = 104;
double thetamin = 0., thetamax = 0.12;

const int nXYBin = 400;
double xymin = -50., xymax = 50.;

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;

AliFOCALGeometry *fGeom;

TClonesArray *fTracks;
TClonesArray *fClusters;
TClonesArray *fTrackClusters = new TClonesArray("AliJBaseTrack", 1500);
TClonesArray *fClustersMatched = new TClonesArray("AliJBaseTrack", 1500);
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
TH2D *hPhiTheta;
TH2D *hXY;
TH2D *hEtaEff;

TH2D *hPhiEtaGamma;
TH2D *hPhiThetaGamma;
TH2D *hXYGamma;

TH1D *hEPhotonCluster;
TH1D *hEPhotonTrue;

Int_t LoadInput(TString inputfile);
void InitOutput();
void FillTruePions();
void FillRecPions(TClonesArray *clusters);
int GetBin(double arr[], int nArr, double val);
AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2);
