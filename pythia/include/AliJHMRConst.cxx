#include "AliJHMRConst.h"

const double detEta[NDET][2] = {
    {3.5, 5.8},
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

const double massWindowMin = 110.;
const double massWindowMax = 160.;
const double massSigmaTrigg[NTRIGGBINS] = {19.3765, 15.2168, 15.0705, 14.2438};
const double massSigmaAssoc[NASSOCBINS] = {21.0777, 29.3914, 19.3765, 15.2168, 15.0705};
const double massPeakPosTrigg[NTRIGGBINS] = {135.685, 140.909, 141.813, 143.411};
const double massPeakPosAssoc[NASSOCBINS] = {139.129, 130., 135.685, 140.909, 141.813};

const double pi0eff = 0.9795;
const double effCorrTrigg[NTRIGGBINS] = {0.221901, 0.409597, 0.440636, 0.42808};
const double effCorrAssoc[NASSOCBINS] = {0.18273, 0.217738, 0.221901, 0.409597, 0.440636};

const double etacut = 0.2;
const double asymcut = 0.7;

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
