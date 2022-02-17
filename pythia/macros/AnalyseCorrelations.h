const int ndata_star = 2;
const int ndata_focal = 1;
const int nTriggBins = 2;
const int nAssocBins = 3;
const double triggPt[nTriggBins+1] = {2.0, 2.5, 3.0};
const double assocPt[nAssocBins+1] = {1.0, 1.5, 2.0, 2.5};
//const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
//const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

// Input
int nEvent;
int nRealTrigg[nTriggBins];

TH1D *hCounter;
TH1D *hRealTriggCounter;
TH1D *hMassTrigg[nTriggBins];
TH1D *hMassAssocPeak[nTriggBins][nAssocBins];
TH1D *hMassAssocSide[nTriggBins][nAssocBins];

TH2D *hCorrReal[nTriggBins][nAssocBins];
TH2D *hCorrMeas[nTriggBins][nAssocBins];
TH2D *hCorrMassMass[nTriggBins][nAssocBins];
TH2D *hCorrMassSide[nTriggBins][nAssocBins];
TH2D *hCorrSideMass[nTriggBins][nAssocBins];
TH2D *hCorrSideSide[nTriggBins][nAssocBins];

// Output
double parTrigg[8];
double parAssoc[9];

double alpha[nTriggBins][nAssocBins];
double beta[nTriggBins];
double st[nTriggBins];

TF1 *fFitTrigg[nTriggBins];
TF1 *fPeakTrigg[nTriggBins];
TF1 *fBgTrigg[nTriggBins];

TF1 *fFitAssoc[nTriggBins][nAssocBins];
TF1 *fPeakAssoc[nTriggBins][nAssocBins];
TF1 *fBgAssoc[nTriggBins][nAssocBins];

TH1D *hCorrRealProj[nTriggBins][nAssocBins];
TH1D *hCorrMeasProj[nTriggBins][nAssocBins];
TH1D *hCorrMassMassProj[nTriggBins][nAssocBins];
TH1D *hCorrMassSideProj[nTriggBins][nAssocBins];
TH1D *hCorrSideMassProj[nTriggBins][nAssocBins];
TH1D *hCorrSideSideProj[nTriggBins][nAssocBins];
TH1D *hCorr[nTriggBins][nAssocBins];

TFile *fIn;
TFile *fOut;

TCanvas *cMassTrigg[nTriggBins];
TCanvas *cMassAssocPeak[nTriggBins][nAssocBins];

// For drawing purposes
TF1 *fPeakColored[nTriggBins];
TF1 *fLeftSidebandColored[nTriggBins];
TF1 *fRightSidebandColored[nTriggBins];
TF1 *fPeakColoredAssoc[nTriggBins][nAssocBins];
TF1 *fLeftSidebandColoredAssoc[nTriggBins][nAssocBins];
TF1 *fRightSidebandColoredAssoc[nTriggBins][nAssocBins];

void processDataSTAR();
void processDataFoCal();

void LoadInput();
void DoAnalysis();
void FitMassPeaks();
void GetScaleFactors();

double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);

void DrawMassHistos(TString dataname);
