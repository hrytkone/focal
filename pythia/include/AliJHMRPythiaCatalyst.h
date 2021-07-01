//===========================================================
// AliJHMRPythiaCatalyst.h
// DongJo Kim (dong.jo.kim@cern.ch)
//===========================================================

#ifndef ALIJHMRPYTHIACATALYST_H
#define ALIJHMRPYTHIACATALYST_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <Pythia8/Pythia.h>

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TRandom3.h"

#include "AliJBaseTrack.h"
#include "AliJHMRHist.h"
#include "AliJHMRConst.h"
#include "set"
#include "map"

using namespace std;
using namespace Pythia8;

class AliJHMRPythiaCatalyst {

	public:

		AliJHMRPythiaCatalyst (Event &inevent, AliJHMRHist *inhistos): 
			event(inevent),
			histos(inhistos){

				unif = new TRandom3();
				fInputList = new TClonesArray("AliJBaseTrack", 1500 );
				TrackEtaRange = 0.8;
			}

		void InitializeEvent(){
			fInputList->Clear();
			UniqueID=0;
		}

		void GetParticles(particleType iType);
		TClonesArray * GetInputList() const {return fInputList;}

		TRandom3 *unif;

		int UniqueID;
		TLorentzVector lvParticle;
		TClonesArray *fInputList;  // tracklist

		Event &event;
		AliJHMRHist *histos;

		double TrackEtaRange ;

};

#endif

