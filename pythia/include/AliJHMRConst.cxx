#include "AliJHMRConst.h"

const double detEta[NDET][2] = {
    {3.5, 5.5},
	{2.6, 4.},
	{-0.8, 0.8},
	{-1000., 1000.}
};

// STAR bins
//double triggPt[NTRIGGBINS+1] = {1.0, 2.0, 2.5, 3.0};
//double assocPt[NASSOCBINS+1] = {0.5, 1.0, 1.5, 2.0, 2.5};

// FoCal bins
double triggPt[NTRIGGBINS+1] = {4.0, 8.0, 20.0};
double assocPt[NASSOCBINS+1] = {2.0, 3.0, 4.0, 8.0};
//double triggPt[NTRIGGBINS+1] = {4.0, 8.0, 10.0, 15.0, 20.0};
//double assocPt[NASSOCBINS+1] = {2.0, 3.0, 4.0, 8.0, 10.0, 15.0};

// When using leading trigger
double leadingPt[NLEADINGBINS+1] = {2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 15., 20.}; // Needed for corrections

// For pure PYTHIA simulations
const double massWindowMin = 110.;
const double massWindowMax = 160.;

// For full detector sim
//const double massWindowMin = 50.;
//const double massWindowMax = 220.;

const double sidebandMin = 300.;
const double sidebandMax = 450.;
const double pi0eff = 0.9795;

const double massSigmaTrigg[NTRIGGBINS] = {27.4508, 28.5018};//, 24.3633, 20.5244};
const double massSigmaAssoc[NASSOCBINS] = {28.8498, 25.1395};//, 27.4508, 28.5018, 24.3633};
const double massPeakPosTrigg[NTRIGGBINS] = {140.779, 143.724};//, 146.13, 146.865};
const double massPeakPosAssoc[NASSOCBINS] = {130, 135.607};//, 140.779, 143.724, 146.13};
const double effCorrAssoc[NASSOCBINS] = {0.297458, 0.383115};//, 0.485133, 0.610642, 0.674529};

// Use these for leading trigger case
const double effCorrLeading[NLEADINGBINS] = {0.250451, 0.275956, 0.326522, 0.368184, 0.424019, 0.527646, 0.546536, 0.585381, 0.618371, 0.620072};

const double etacut = 0.2;
const double asymcut = 1.0;
const double thetacut = 0.;

const double kJPi = TMath::Pi();
const double kJTwoPi = 2 * TMath::Pi();

const double etaBinWidth = 0.025;
const double phiBinWidth = 0.025;

const double etaTrackerRange = 0.9;
const double etaFocalRange = detEta[0][1] - detEta[0][0];
const double etaSTARRange = detEta[1][1] - detEta[1][0];
const int nEtaBinTracker = 2.*int(etaTrackerRange/etaBinWidth) + 1;
const int nEtaBinFocal = int(etaFocalRange/etaBinWidth) + 1;
const int nEtaBinSTAR = int(etaSTARRange/etaBinWidth) + 1;

const double deltaPhiMin = -kJPi/2.0;
const double deltaPhiMax = 3.0*kJPi/2.0;
const int nPhiBin = int((deltaPhiMax-deltaPhiMin)/phiBinWidth) + 1;

double logBinsX[NINCPTBIN+1], limMin = 0.1, limMax = 100;
const double logBW = (log(limMax) - log(limMin))/NINCPTBIN;

const double incEtaRange = 20.0;
double limPhotonEnergyMin = 0., limPhotonEnergyMax = 1500.;
