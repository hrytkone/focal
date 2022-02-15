#ifndef ALIJHMRGEANTCATALYST_H
#define ALIJHMRGEANTCATALYST_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"

#include "AliJBaseTrack.h"
#include "AliJHMREvent.h"
#include "AliJHMRCluster.h"
#include "AliJHMRTrack.h"
#include "AliJHMRHist.h"
#include "AliJHMRConst.h"
#include "set"
#include "map"

using namespace std;

class AliJHMRGeantCatalyst {

	public:

		AliJHMRGeantCatalyst (TString infile, AliJHMRHist *inhistos) :
			inputfile(infile),
			histos(inhistos) {
                fInputListCluster = new TClonesArray("AliJBaseTrack", 1500);
				fInputListTrack   = new TClonesArray("AliJBaseTrack", 1500);
                fEvent = 0;
                fClusters = 0;
                fTracks = 0;
		}

        Int_t LoadInput();
        Int_t GetNumberOfEvents() { return fTree->GetEntries(); }
        void GetEvent(int iev) { fTree->GetEntry(iev); }
        TClonesArray * GetClusters();
        TClonesArray * GetPi0True();

        int UniqueID;

        TLorentzVector lvParticle;
        TString inputfile;
        TFile *fIn;
        TTree *fTree;
        TBranch *fBrtrack;
        TBranch *fBrcluster;
        AliJHMREvent *fEvent;

		AliJHMRHist *histos;
        TClonesArray *fClusters;
        TClonesArray *fTracks;
        TClonesArray *fInputListCluster;
        TClonesArray *fInputListTrack;
};

#endif
