const TString outputname = "accceptance_check";
const double asymcut = 1.;

const int nIncPtBin = 150;
double logBinsX[nIncPtBin+1], limMin = 0.1, limMax = 100;
double logBW = (log(limMax) - log(limMin))/nIncPtBin;

const double etamin = 2.8;
const double etamax = 6.2;
const double etaRange = etamax - etamin;
const double etaBinWidth = 0.025;
int nEtaBin = int(etaRange/etaBinWidth) + 1;

const double xmin = -60.;
const double xmax = 60.;
const double xRange = xmax - xmin;
const double xBinWidth = 0.025;
int nXBin = int(xRange/xBinWidth) + 1;

const double xcutmin = -3.;
const double xcutmax = 3.;
const double ycutmin = -3.;
const double ycutmax = 3.;

TFile *fIn, *fOut;
TTree *fTree;
TBranch *fBrtrack;
TBranch *fBrcluster;
AliJHMREvent *fEvent;

AliFOCALGeometry *fGeom;

TClonesArray *fTracks;
TClonesArray *fClusters;
TLorentzVector lvParticle;

// Histograms
TH2D *hXEta;
TH2D *hYEta;
TH2D *hXY;

Int_t LoadInput(TString inputfile);
void InitOutput();
