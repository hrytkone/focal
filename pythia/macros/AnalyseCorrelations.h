const bool bUseConstMassWindow = true;
const bool useLeading = false;

const int ndata_star = 2;
const int ndata_focal = 1;
const int nTriggBins = 2;
const int nAssocBins = 2;
const int nLeadingBins = 10;

const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};
//const double triggPt[nTriggBins+1] = {8.0, 20.0};
//const double assocPt[nAssocBins+1] = {3.0, 4.0, 8.0};
//const double triggPt[nTriggBins+1] = {4.0, 8.0, 10.0, 15.0, 20.0};
//const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0, 8.0, 10.0, 15.0};

// Leading trigger
//const double triggPt[nTriggBins+1] = {4.0, 10000.0};
//const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

const double pi0br = 1./0.98823;
//const double pi0eff = 1.;

const double effCorrTrigg[nTriggBins] = {1., 1.};
//const double effCorrTrigg[nTriggBins] = {1.};
//const double effCorrTrigg[nTriggBins] = {1./0.485133, 1./0.610642, 1./0.674529, 1./0.561265};
//const double effCorrTrigg[nTriggBins] = {1./0.221901, 1./0.409597, 1./0.440636, 1./0.42808};
//const double effCorrTrigg[nTriggBins] = {1./0.988, 1./0.988};
//const double effCorrTrigg[nTriggBins] = {1./0.43981}; // LEADING TRIGG 2 GeV/c < pTassoc < pTtrigg

const double massMin = 110.;
const double massMax = 160.;

//const double massMin = 60.;
//const double massMax = 210.;
const double massSigmaTrigg[nTriggBins] = {27.4508};//, 28.5018, 24.3633, 20.5244};
const double massSigmaAssoc[nAssocBins] = {28.8498, 25.1395};//, 27.4508, 28.5018, 24.3633};
const double massPeakPosTrigg[nTriggBins] = {140.779};//, 143.724, 146.13, 146.865};
const double massPeakPosAssoc[nAssocBins] = {130, 135.607};//, 140.779, 143.724, 146.13};

// Input
int nEvent;
int nRealTrigg[nTriggBins];

TH1D *hCounter;
TH1D *hRealTriggCounter;
TH1D *hMassTrigg[nTriggBins];
TH1D *hMassAssocPeak[nTriggBins][nAssocBins];
TH1D *hMassAssocSide[nTriggBins][nAssocBins];
TH1D *hMassAssocSum[nTriggBins][nAssocBins];

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
