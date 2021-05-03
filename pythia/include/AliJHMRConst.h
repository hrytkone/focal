    //-------------------------
    //         pT bins        |
    //-------------------------
    const int nTriggBins = 4;
    double  triggPt[nTriggBins+1] = {1.0, 2.0, 4.0, 8.0, 20.0};
    //double  triggPt[nTriggBins+1] = {3.0, 1000.0};
    //double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    const int nAssocBins = 4;
    double  assocPt[nAssocBins+1] = {0.5, 1.0, 2.0, 3.0, 4.0};
    //double  assocPt[nAssocBins+1] = {3.0, 1000.0};
    //double  assocPt[nAssocBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};


    const double etaBinWidth = 0.025;
    const double phiBinWidth = 0.025;

    const double etaTrackerRange = 0.9;
    const double etaFocalMin = 3.2;
    const double etaFocalMax = 5.8;
    const double etaFocalRange = etaFocalMax - etaFocalMin;
    const int nEtaBinTracker = 2.*int(etaTrackerRange/etaBinWidth) + 1;
    const int nEtaBinFocal = int(etaFocalRange/etaBinWidth) + 1;

    const double deltaPhiMin = -TMath::Pi()/2.0;
    const double deltaPhiMax = 3.0*TMath::Pi()/2.0;
    const int nPhiBin = int((deltaPhiMax-deltaPhiMin)/phiBinWidth) + 1;

    const int nIncPtBin = 150;
    double logBinsX[nIncPtBin+1], limMin = 0.1, limMax = 100;
    const double logBW = (log(limMax) - log(limMin))/nIncPtBin;

    const int nIncEtaBin = 150;
    const double incEtaRange = 20.0;

    const int nPhotonEnergyBin = 150;
    double limPhotonEnergyMin = 0., limPhotonEnergyMax = 1500.;