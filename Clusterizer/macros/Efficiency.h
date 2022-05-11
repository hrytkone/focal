const TString outputname = "efficiency_small-bins_no-asym-cut";
const double asymcut = 1.;

const int nIncPtBin = 150;
double logBinsX[nIncPtBin+1], limMin = 0.1, limMax = 100;
double logBW = (log(limMax) - log(limMin))/nIncPtBin;

const double etamin = 3.2;
const double etamax = 5.8;
const double etaRange = etamax - etamin;
const double etaBinWidth = 0.025;
int nEtaBin = int(etaRange/etaBinWidth) + 1;

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

TH1D *hEPhotonCluster;
TH1D *hEPhotonTrue;

Int_t LoadInput(TString inputfile);
void InitOutput();
void FillTruePions();
void FillRecPions();
AliJBaseTrack GetPhotonSumVector(TClonesArray *clusters, AliJBaseTrack *lv1, AliJBaseTrack *lv2);
