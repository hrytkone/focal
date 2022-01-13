const int ndata = 2;
const int nTriggBins = 2;
const int nAssocBins = 3;
const double triggPt[nTriggBins+1] = {2.0, 2.5, 3.0};
const double assocPt[nAssocBins+1] = {1.0, 1.5, 2.0, 2.5};

// Input
int nEvent;
int nRealTrigg[nTriggBins];

TH1D *hCounter;
TH1D *hRealTriggCounter;
TH1D *hMassTrigg[nTriggBins];
TH1D *hMassAssocPeak[nTriggBins][nAssocBins];
TH1D *hMassAssocSide[nTriggBins][nAssocBins];

TH2D *hCorrReal[nTriggBins][nAssocBins];
TH2D *hCorrMassMass[nTriggBins][nAssocBins];
TH2D *hCorrMassSide[nTriggBins][nAssocBins];
TH2D *hCorrSideMass[nTriggBins][nAssocBins];
TH2D *hCorrSideSide[nTriggBins][nAssocBins];

// Output
double alphat[nTriggBins];
double betat[nTriggBins];
double alphaa[nTriggBins][nAssocBins];
double betaa[nTriggBins][nAssocBins];
double A[nTriggBins][nAssocBins], B[nTriggBins][nAssocBins];

TF1 *fFitTrigg[nTriggBins];
TF1 *fPeakTrigg[nTriggBins];
TF1 *fBgTrigg[nTriggBins];
TF1 *fFitAssoc[nTriggBins][nAssocBins];
TF1 *fPeakAssoc[nTriggBins][nAssocBins];
TF1 *fBgAssoc[nTriggBins][nAssocBins];

TH1D *hCorrRealProj[nTriggBins][nAssocBins];
TH1D *hCorrMassMassProj[nTriggBins][nAssocBins];
TH1D *hCorrMassSideProj[nTriggBins][nAssocBins];
TH1D *hCorrSideMassProj[nTriggBins][nAssocBins];
TH1D *hCorrSideSideProj[nTriggBins][nAssocBins];
TH1D *hCorr[nTriggBins][nAssocBins];

TFile *fIn;
TFile *fOut;

TCanvas *cMassesTrigg[ndata];
TCanvas *cMassesAssoc[ndata];

void processDataSTAR();

void LoadInput();
void DoAnalysis();
void GetNormalisations();

double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);

void PlotMassHistos();
