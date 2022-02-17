#ifndef ALIJHMRCONST_H
#define ALIJHMRCONST_H

#include "TMath.h"
#include "TString.h"

#define NDET 4
#define NTRIGGBINS 3
#define NASSOCBINS 4
#define NINCPTBIN 150
#define NINCETABIN 150
#define NPHOTONENERGYBIN 150

enum particleType { kJHadron, kJPi0, kJDecayPhoton, kJRecPi0 };
enum detector { kJFoCal, kJSTAR, kJTracker, kJFull };

extern const double etacut;

extern const double detEta[NDET][2];
extern const TString accFunc[NDET];

extern double triggPt[NTRIGGBINS+1];
extern double assocPt[NASSOCBINS+1];

extern const double kJPi;
extern const double kJTwoPi;

extern const double etaBinWidth;
extern const double phiBinWidth;

extern const double etaTrackerRange;
extern const double etaFocalRange;
extern const double etaSTARRange;
extern const int nEtaBinTracker;
extern const int nEtaBinFocal;
extern const int nEtaBinSTAR;

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
