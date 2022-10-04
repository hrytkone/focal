#include "AliJHMRConst.h"

const double detEta[NDET][2] = {
    {4.2, 5.3},
    //{3.5, 5.8},
	{2.6, 4.},
	{-0.8, 0.8},
	{-1000., 1000.}
};

const TString accFunc[NDET] = {
    "TMath::Exp(-0.117082/(x + 0.0832931))", // UPDATE THIS
    "TMath::Exp(-0.123048/(x + 0.0874463))",
    "TMath::Exp(-0.117082/(x + 0.0832931))", // NOT CORRECT, CHECK CORRECT ONE
    "1."
};

// STAR bins
//double triggPt[NTRIGGBINS+1] = {1.0, 2.0, 2.5, 3.0};
//double assocPt[NASSOCBINS+1] = {0.5, 1.0, 1.5, 2.0, 2.5};

// FoCal bins
//double triggPt[NTRIGGBINS+1] = {4.0, 8.0, 20.0};
//double assocPt[NASSOCBINS+1] = {2.0, 3.0, 4.0};
double triggPt[NTRIGGBINS+1] = {4.0, 8.0, 10.0, 15.0, 20.0};
double assocPt[NASSOCBINS+1] = {2.0, 3.0, 4.0, 8.0, 10.0, 15.0};

// For pure PYTHIA simulations
//const double massWindowMin = 110.;
//const double massWindowMax = 160.;

const double massWindowMin = 50.;
const double massWindowMax = 220.;

const double sidebandMin = 300.;
const double sidebandMax = 450.;

const double massSigmaTrigg[NTRIGGBINS] = {16.896, 16.8545, 17.8731, 11.8696};
const double massSigmaAssoc[NASSOCBINS] = {26.3631, 20.0558, 16.896, 16.8545, 17.8731};
const double massPeakPosTrigg[NTRIGGBINS] = {140.115, 142.753, 144.739, 143.427};
const double massPeakPosAssoc[NASSOCBINS] = {130.582, 135.952, 140.115, 142.753, 144.739};

const double pi0eff = 0.9795;
//const double effCorrTrigg[NTRIGGBINS] = {0.47261, 0.590262, 0.649246, 0.477833};
//const double effCorrAssoc[NASSOCBINS] = {0.2903, 0.385151, 0.47261, 0.590262, 0.649246};
const double effCorrTrigg[NTRIGGBINS] = {0.485133, 0.610642, 0.674529, 0.561265};
const double effCorrAssoc[NASSOCBINS] = {0.297458, 0.383115, 0.485133, 0.610642, 0.674529};


const double etacut = 0.;
const double asymcut = 0.8;

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
