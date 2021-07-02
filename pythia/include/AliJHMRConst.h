#ifndef ALIJHMRCONST_H
#define ALIJHMRCONST_H

#include "TMath.h"

#define NTRIGGBINS 4
#define NASSOCBINS 4
#define NINCPTBIN 150
#define NINCETABIN 150
#define NPHOTONENERGYBIN 150

enum particleType { kJHadron, kJPi0, kJDecayPhoton };

//extern const int nTriggBins;
extern double triggPt[NTRIGGBINS+1];

//extern const int nAssocBins;
extern double assocPt[NASSOCBINS+1];

extern const double kJPi;
extern const double kJTwoPi;

extern const double etaBinWidth;
extern const double phiBinWidth;

extern const double etaTrackerRange;
extern const double etaFocalMin;
extern const double etaFocalMax;
extern const double etaFocalRange;
extern const int nEtaBinTracker;
extern const int nEtaBinFocal;

extern const double deltaPhiMin;
extern const double deltaPhiMax;
extern const int nPhiBin;

//extern const int nIncPtBin;
extern double logBinsX[NINCPTBIN+1], limMin, limMax;
extern const double logBW;

//extern const int nIncEtaBin = 150;
extern const double incEtaRange;

//extern const int nPhotonEnergyBin = 150;
extern double limPhotonEnergyMin, limPhotonEnergyMax;

#endif