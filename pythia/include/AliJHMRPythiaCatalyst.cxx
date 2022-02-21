#include "AliJHMRPythiaCatalyst.h"
#include <set>

using namespace std;

void AliJHMRPythiaCatalyst::GetParticles(detector idet) {

	for (int partIdx = 0; partIdx < event.size(); partIdx++) {

		double trEta = event[partIdx].eta();
		if( trEta < detEta[idet][0] || trEta > detEta[idet][1] ) continue;

		lvParticle.SetPxPyPzE(event[partIdx].px(), event[partIdx].py(), event[partIdx].pz(), event[partIdx].e() );
		lvParticle.SetUniqueID(UniqueID++);

        Int_t motherId = event[partIdx].mother1();

		AliJBaseTrack track( lvParticle );
        track.SetID(event[partIdx].id());
		track.SetMCIndex(partIdx);
        track.SetMotherID(motherId);
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

            // In the case of photons tag those that are from pi0 decay
            //      track label 1 = decay product
            //      track label 0 = not from decay
            if (event[motherId].id() == 111 && (event[motherId].eta() > detEta[idet][0] && event[motherId].eta() < detEta[idet][1])) {
                track.SetLabel(1);
            } else {
                track.SetLabel(0);
            }
			new((*fInputListPhoton)[fInputListPhoton->GetEntriesFast()]) AliJBaseTrack(track);
		}
	}
}

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
