const TString outputname = "efficiency_pythiamb_asym-08";
const double asymcut = 0.8;

const int nPtBin = 6;
double pt[nPtBin+1], limMin = 2, limMax = 20;
double logBW = (log(limMax) - log(limMin))/nPtBin;

const int nEtaBin = 56;
//const int nEtaBin = 30;
double eta[nEtaBin+1];
double etaBW = 0.05, etamin = 3.0, etamax = 5.8;
//double etaBW = 0.05, etamin = 3.2, etamax = 5.5;

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
