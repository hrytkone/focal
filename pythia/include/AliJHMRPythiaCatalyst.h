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

		AliJHMRPythiaCatalyst (Event &inevent, AliJHMRHist *inhistos) : 
			event(inevent),
			histos(inhistos){
				unif = new TRandom3();
				fInputListHadron = new TClonesArray("AliJBaseTrack", 1500);
				fInputListPi0 = new TClonesArray("AliJBaseTrack", 1500);
				fInputListPhoton = new TClonesArray("AliJBaseTrack", 1500);								
			}

		void InitializeEvent(){
			fInputListHadron->Clear("C");			
			fInputListPi0->Clear("C");			
			fInputListPhoton->Clear("C");
			UniqueID=0;
		}

		void GetParticles(detector idet);
		TClonesArray * GetParticleList(particleType itype);

		TRandom3 *unif;

		int UniqueID;
		TLorentzVector lvParticle;
		TClonesArray *fInputListHadron;
		TClonesArray *fInputListPi0;	
		TClonesArray *fInputListPhoton;

		Event &event;
		AliJHMRHist *histos;
};

#endif

