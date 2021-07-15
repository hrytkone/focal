#include "AliJHMRConst.h"

const double detEta[NDET][2] = {
	{3.2, 5.8},
	{-0.8, 0.8},
	{-1000., 1000.}
};

double triggPt[NTRIGGBINS+1] = {1.0, 2.0, 4.0, 8.0, 20.0};
double assocPt[NASSOCBINS+1] = {0.5, 1.0, 2.0, 3.0, 4.0};

const double kJPi = TMath::Pi();
const double kJTwoPi = 2 * TMath::Pi();

const double etaBinWidth = 0.025;
const double phiBinWidth = 0.025;

const double etaTrackerRange = 0.9;
const double etaFocalMin = 3.2;
const double etaFocalMax = 5.8;
const double etaFocalRange = etaFocalMax - etaFocalMin;
const int nEtaBinTracker = 2.*int(etaTrackerRange/etaBinWidth) + 1;
const int nEtaBinFocal = int(etaFocalRange/etaBinWidth) + 1;

const double deltaPhiMin = -kJPi/2.0;
const double deltaPhiMax = 3.0*kJPi/2.0;
const int nPhiBin = int((deltaPhiMax-deltaPhiMin)/phiBinWidth) + 1;

double logBinsX[NINCPTBIN+1], limMin = 0.1, limMax = 100;
const double logBW = (log(limMax) - log(limMin))/NINCPTBIN;

const double incEtaRange = 20.0;
double limPhotonEnergyMin = 0., limPhotonEnergyMax = 1500.;