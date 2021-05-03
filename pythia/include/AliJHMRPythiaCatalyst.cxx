#include "AliJHMRPythiaCatalyst.h"
#include <set>

using namespace std;

//extern double DeltaPhiX(double dphi);

//--------------------------------------------------------
void AliJHMRPythiaCatalyst::GetParticles(particleType iType) { 

	if (iType==kJHadron)

		for(int partIdx = 0; partIdx < event.size(); partIdx++){ 

			if( !event[partIdx].isFinal() || !event[partIdx].isCharged() || !event[partIdx].isHadron() ) continue;
			//if(event[partIdx].vProd().pAbs()>1.0) continue; //production vertex >1.0mm to reject secondaries. 

			double trEta = event[partIdx].eta();
			//cout<< "pTt=" << pTtBin <<" pTa="<<pTaBin <<endl;  
			if( fabs(trEta)>TrackEtaRange ) continue;

			lvParticle.SetPxPyPzE(event[partIdx].px(), event[partIdx].py(), event[partIdx].pz(), event[partIdx].e() );
			lvParticle.SetUniqueID(UniqueID++);
			// Make a jbasetrack for particle arrau
			AliJBaseTrack track( lvParticle );
			track.SetID(event[partIdx].id());
			track.SetParticleType(kJHadron);
			track.SetCharge(event[partIdx].charge());
			track.SetTrackEff(1.);

			new((*fInputList)[fInputList->GetEntriesFast()]) AliJBaseTrack(track);
		} // end of particle loop

	}

} // end of function

