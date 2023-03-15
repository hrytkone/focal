const int ndata_star = 2;
const int ndata_focal = 1;
const int nTriggBins = 2;
const int nAssocBins = 2;
//const int nTriggBins = 3;
//const int nAssocBins = 4;

const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};
//const double triggPt[nTriggBins+1] = {1.0, 2.0, 2.5, 3.0};
//const double assocPt[nAssocBins+1] = {0.5, 1.0, 1.5, 2.0, 2.5};

const double pi0br = 1./0.98823;
//const double pi0eff = 1.;

const bool bPythia = true;
const double massMin = 120.;
const double massMax = 155.;
//const double massMin = 50.;
//const double massMax = 200.;

const double sbMin = 250.;
const double sbMax = 450.;

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
TH2D *hCorrMeasMixed[nTriggBins][nAssocBins];

// Same event
TH2D *hCorrMassMass[nTriggBins][nAssocBins];
TH2D *hCorrMassSide[nTriggBins][nAssocBins];
TH2D *hCorrSideMass[nTriggBins][nAssocBins];
TH2D *hCorrSideSide[nTriggBins][nAssocBins];

// Mixed event
TH2D *hCorrMassMassMixed[nTriggBins][nAssocBins];
TH2D *hCorrMassSideMixed[nTriggBins][nAssocBins];
TH2D *hCorrSideMassMixed[nTriggBins][nAssocBins];
TH2D *hCorrSideSideMixed[nTriggBins][nAssocBins];

// Output
double parTrigg[8];
double parAssocPeak[9];
double parAssocSide[9];
double parAssocSum[9];
double fitConstantVal[nTriggBins] = {0};

double triggMean[nTriggBins] = {0};
double assocMeanPeak[nTriggBins][nAssocBins] = {0};
double assocMeanSide[nTriggBins][nAssocBins] = {0};
double triggSigma[nTriggBins] = {0};
double assocSigmaPeak[nTriggBins][nAssocBins] = {0};
double assocSigmaSide[nTriggBins][nAssocBins] = {0};

double alpha[nTriggBins][nAssocBins];
double beeta[nTriggBins][nAssocBins];
double yamma[nTriggBins][nAssocBins];
double st[nTriggBins];
double stobtrigg[nTriggBins];
double stobassoc[nTriggBins][nAssocBins];
double stfit[nTriggBins];

TF1 *fFitTrigg[nTriggBins];
TF1 *fPeakTrigg[nTriggBins];
TF1 *fBgTrigg[nTriggBins];

TF1 *fFitAssocPeak[nTriggBins][nAssocBins];
TF1 *fPeakAssocPeak[nTriggBins][nAssocBins];
TF1 *fBgAssocPeak[nTriggBins][nAssocBins];

TF1 *fFitAssocSide[nTriggBins][nAssocBins];
TF1 *fPeakAssocSide[nTriggBins][nAssocBins];
TF1 *fBgAssocSide[nTriggBins][nAssocBins];

TF1 *fFitAssocSum[nTriggBins][nAssocBins];
TF1 *fPeakAssocSum[nTriggBins][nAssocBins];
TF1 *fBgAssocSum[nTriggBins][nAssocBins];

TH1D *hCorrRealProj[nTriggBins][nAssocBins];
TH1D *hCorrMassMassProj[nTriggBins][nAssocBins];
TH1D *hCorrMassSideProj[nTriggBins][nAssocBins];
TH1D *hCorrSideMassProj[nTriggBins][nAssocBins];
TH1D *hCorrSideSideProj[nTriggBins][nAssocBins];
TH1D *hCorrNonCorrected[nTriggBins][nAssocBins];
TH1D *hCorrCorrected[nTriggBins][nAssocBins];

TFile *fIn;
TFile *fOut;

TCanvas *cMassTrigg[nTriggBins];
TCanvas *cMassAssocPeak[nTriggBins][nAssocBins];
TCanvas *cMassAssocSide[nTriggBins][nAssocBins];

// For drawing purposes
TF1 *fPeakColored[nTriggBins];
TF1 *fBgColored[nTriggBins];
TF1 *fLeftSidebandColored[nTriggBins];
TF1 *fRightSidebandColored[nTriggBins];
TF1 *fPeakColoredAssocPeak[nTriggBins][nAssocBins];
TF1 *fLeftSidebandColoredAssocPeak[nTriggBins][nAssocBins];
TF1 *fRightSidebandColoredAssocPeak[nTriggBins][nAssocBins];
TF1 *fBgColoredAssocPeak[nTriggBins][nAssocBins];
TF1 *fBgColoredAssocSide[nTriggBins][nAssocBins];
TF1 *fPeakColoredAssocSide[nTriggBins][nAssocBins];
TF1 *fLeftSidebandColoredAssocSide[nTriggBins][nAssocBins];
TF1 *fRightSidebandColoredAssocSide[nTriggBins][nAssocBins];

void processDataSTAR();
void processDataFoCal();

void LoadInput();
void DoAnalysis();
void DoSimplifiedAnalysis();
void MixedEventCorrection();
double GetMixedEventNormalization(TH2D *h);
void FitMassPeaks();
void GetScaleFactorsVersion1();
void GetScaleFactorsVersion2();
void GetScaleFactorsVersion3();

double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);

void DrawMassHistos(TString dataname);
void SetStyle(Bool_t graypalette);
void redrawBorder();
