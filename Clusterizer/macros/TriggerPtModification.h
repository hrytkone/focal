const TString outputname = "triggPt-mod_binning";
const double asymcut = 0.8;

const double etaMin = 3.5;
const double etaMax = 5.5;
const double etacut = 0.2;

const double mwMin = 110;
const double mwMax = 160;

const int nPtBin = 90;
const double ptMin = 2.;
const double ptMax = 20.;

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;

AliFOCALGeometry *fGeom;

TClonesArray *fTracks;
TClonesArray *fClusters;
TClonesArray *fTruePi0 = new TClonesArray("AliJBaseTrack", 1500);
TClonesArray *fRecPi0 = new TClonesArray("AliJBaseTrack", 1500);
TLorentzVector lvParticle;

// Histograms
TH1D *hPtTrue;
TH1D *hPtRec;

Int_t LoadInput(TString inputfile);
void InitOutput();
void GetTruePions();
void GetRecPions();
void FillTruePt();
void FillRecPt();
AliJBaseTrack GetPhotonSumVector(AliJBaseTrack *lv1, AliJBaseTrack *lv2);
double GetAsymmetry(double E1, double E2);
