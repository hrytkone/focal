#include "AliJHMRPythiaCatalyst.h"
#include <set>

using namespace std;

void AliJHMRPythiaCatalyst::GetParticles(detector idet) { 

		for (int partIdx = 0; partIdx < event.size(); partIdx++) { 

			double trEta = event[partIdx].eta();
			if( trEta < detEta[idet][0] || trEta > detEta[idet][1] ) continue;

			lvParticle.SetPxPyPzE(event[partIdx].px(), event[partIdx].py(), event[partIdx].pz(), event[partIdx].e() );
			lvParticle.SetUniqueID(UniqueID++);
			
			AliJBaseTrack track( lvParticle );
			track.SetID(event[partIdx].id());
			track.SetCharge(event[partIdx].charge());
			track.SetTrackEff(1.);

			if ( event[partIdx].isFinal() && event[partIdx].isCharged() && event[partIdx].isHadron() ) {
				track.SetParticleType(kJHadron);	
				new((*fInputListHadron)[fInputListHadron->GetEntriesFast()]) AliJBaseTrack(track);
			}

			if ( event[partIdx].id() == 111 ) {
				track.SetParticleType(kJPi0);	
				new((*fInputListPi0)[fInputListPi0->GetEntriesFast()]) AliJBaseTrack(track);
			}

			if ( event[partIdx].isFinal() && event[partIdx].id() == 22 ) {
				track.SetParticleType(kJDecayPhoton);	
				new((*fInputListPhoton)[fInputListPhoton->GetEntriesFast()]) AliJBaseTrack(track);
			}
		} // end of particle loop

} // end of function

TClonesArray * AliJHMRPythiaCatalyst::GetParticleList(particleType itype)
{
	switch(itype) {
		case kJHadron:
			return fInputListHadron;
		case kJPi0:
			return fInputListPi0;
		case kJDecayPhoton:
			return fInputListPhoton;
		default:
			std::cout << "Particle species not specified, return nothing" << std::endl;
			return 0;
	}
}
