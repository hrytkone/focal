const TString outputname = "efficiency_pi0-gun_asym-08";
//const TString outputname = "efficiency_pythiamb";
const double asymcut = 0.8;

//const int nPtBin = 38;
const int nPtBin = 19;
double pt[nPtBin+1], limMin = 2, limMax = 18;
double logBW = (log(limMax) - log(limMin))/nPtBin;

//const int nEtaBin = 38;
const int nEtaBin = 19;
double eta[nEtaBin+1], etamin = 3.4, etamax = 5.3;
//double eta[nEtaBin+1], etamin = 3.0, etamax = 6.0;
double etaBW = (etamax - etamin)/nEtaBin;

//const int nPhiBin = 52;
const int nPhiBin = 104;
double phimin = -TMath::Pi(), phimax = TMath::Pi();

//const int nThetaBin = 52;
const int nThetaBin = 104;
double thetamin = 0., thetamax = 0.12;

const int nXYBin = 400;
double xymin = -50., xymax = 50.;

int ntrue = 0;

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
TH2D *hEtaMass;
TH2D *hPhiEtaTrue;
TH2D *hPhiEta;
TH2D *hPhiTheta;
TH2D *hXY;
TH2D *hEtaEff;

TH2D *hEtaPtRec_match;
TH2D *hEtaERec_match;
TH2D *hPtMass_match;
TH2D *hEtaMass_match;
TH2D *hPhiEta_match;
TH2D *hPhiTheta_match;
TH2D *hXY_match;

TH2D *hEtaPtRec_match_truept;

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
