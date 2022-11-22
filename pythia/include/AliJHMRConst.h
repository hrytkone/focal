#ifndef ALIJHMRCONST_H
#define ALIJHMRCONST_H

#include "TMath.h"
#include "TString.h"

#define NDET 4
#define NTRIGGBINS 4
#define NASSOCBINS 5
//#define NTRIGGBINS 1
//#define NASSOCBINS 1
#define NLEADINGBINS 4
#define NMASSBINS 6
#define NINCPTBIN 150
#define NINCETABIN 150
#define NPHOTONENERGYBIN 150

enum particleType { kJHadron, kJPi0, kJDecayPhoton, kJRecPi0 };
enum detector { kJFoCal, kJSTAR, kJTracker, kJFull };

// For fixed mass window
extern const double massWindowMin;
extern const double massWindowMax;
extern const double sidebandMin;
extern const double sidebandMax;

// mass peak sigmas & peak positions for each pT bin
// DEPENDS ON ETA AND ASYM CUTS
extern const double massSigmaTrigg[NTRIGGBINS];
extern const double massSigmaAssoc[NASSOCBINS];
extern const double massPeakPosTrigg[NTRIGGBINS];
extern const double massPeakPosAssoc[NASSOCBINS];

extern const double pi0eff; // FOR MC SIMULATION
extern const double effCorrTrigg[NTRIGGBINS]; // FOR GEANT SIMUALTION
extern const double effCorrAssoc[NASSOCBINS]; // FOR GEANT SIMUALTION
extern const double effCorrLeadingTrigg[NLEADINGBINS]; // FOR GEANT SIMUALTION

extern const double etacut;
extern const double asymcut;

extern const double detEta[NDET][2];
extern const TString accFunc[NDET];

extern double triggPt[NTRIGGBINS+1];
extern double assocPt[NASSOCBINS+1];
extern double leadingTriggPt[NLEADINGBINS+1];

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
