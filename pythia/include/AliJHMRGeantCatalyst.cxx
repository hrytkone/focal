#include "AliJHMRGeantCatalyst.h"
#include <set>

using namespace std;

void AliJHMRGeantCatalyst::GetClusters()
{
    for (int iclust = 0; iclust < fClusters->GetEntries(); iclust++) {
        AliJHMRCluster *cluster = (AliJHMRCluster*)fClusters->At(iclust);
        //lvParticle.SetPtEtaPhiE(cluster->GetPt(), cluster->GetEta(), cluster->GetPhi(), cluster->GetE());
        lvParticle.SetPtEtaPhiM(cluster->GetPt(), cluster->GetEta(), cluster->GetPhi(), 0.);
        lvParticle.SetUniqueID(UniqueID++);

        AliJBaseTrack track( lvParticle );
        new((*fInputListCluster)[fInputListCluster->GetEntriesFast()]) AliJBaseTrack(track);
    }
    //return fInputListCluster;
}

// Only for foCal
void AliJHMRGeantCatalyst::GetParticles()
{
    double etaMin = detEta[0][0];
    double etaMax = detEta[0][1];
    for (int itrack = 0; itrack < fTracks->GetEntries(); itrack++) {
        AliJBaseTrack *tr = (AliJBaseTrack*)fTracks->At(itrack);
        double trEta = tr->Eta();

        if ( tr->IsPrimary() && tr->GetID() == 22 ) {

            AliJBaseTrack track = *tr;
            track.SetParticleType(kJDecayPhoton);

            // In the case of photons tag those that are from pi0 decay
            //      track label 1 = decay product
            //      track label 0 = not from decay
            Int_t imother = track.GetMotherID();
            AliJBaseTrack *trMom = (AliJBaseTrack*)fTracks->At(imother);
            if (trMom->GetID() == 111 && (trMom->Eta() > etaMin+etacut && trMom->Eta() < etaMax-etacut)) {
                track.SetLabel(1);
            } else {
                track.SetLabel(0); // For now just set to 0
            }
            new((*fInputListPhoton)[fInputListPhoton->GetEntriesFast()]) AliJBaseTrack(track);
        }

        // Use smaller acceptance than for gammas to suppress the effect from missing gamma pairs
        if ( trEta < etaMin+etacut || trEta > etaMax-etacut ) continue;

        if ( tr->GetID() == 111 ) {
            AliJBaseTrack track = *tr;
            track.SetParticleType(kJPi0);
            new((*fInputListPi0)[fInputListPi0->GetEntriesFast()]) AliJBaseTrack(track);
        }
    }
}

void AliJHMRGeantCatalyst::GetPtMatchedClusters(TClonesArray *arrPi0Candidates)
{
    std::vector<bool> bParticleIsMatched;
    for (int i = 0; i < arrPi0Candidates->GetEntriesFast(); i++) bParticleIsMatched.push_back(false);
    for (int irec = 0; irec < arrPi0Candidates->GetEntriesFast(); irec++) {
        AliJBaseTrack *pi0rec = (AliJBaseTrack*)arrPi0Candidates->At(irec);
        double etaRec = pi0rec->Eta();
        double phiRec = pi0rec->Phi();

        // Get geometrically closest MC pi0 to the reconstructed one
        double dist = 1000000.0;
        double ptMatched = -1;
        int iMatched = 0;
        for (int itr = 0; itr < fInputListPi0->GetEntriesFast(); itr++) {
            AliJBaseTrack *tr = (AliJBaseTrack*)fInputListPi0->At(itr);

            double etaTrue = tr->Eta();
            double phiTrue = tr->Phi();
            double ptTrue = tr->Pt();
            double delta = TMath::Sqrt((etaTrue-etaRec)*(etaTrue-etaRec) + (phiTrue-phiRec)*(phiTrue-phiRec));

            if (delta < dist) {
                dist = delta;
                ptMatched = ptTrue;
                if (!bParticleIsMatched[iMatched]) iMatched = itr;
            }
        }
        bParticleIsMatched[iMatched] = true;

        // Save cluster to the list but with matched pt
        lvParticle.SetPtEtaPhiM(ptMatched, pi0rec->Eta(), pi0rec->Phi(), 0.);
        lvParticle.SetUniqueID(UniqueID++);

        AliJBaseTrack track( lvParticle );
        new((*fInputListPtMatchedPi0)[fInputListPtMatchedPi0->GetEntriesFast()]) AliJBaseTrack(track);        
    }
}

TClonesArray * AliJHMRGeantCatalyst::GetParticleList(particleType itype)
{
	switch(itype) {
		case kJHadron:
			return fInputListHadron;
		case kJPi0:
			return fInputListPi0;
		case kJDecayPhoton:
			return fInputListPhoton;
        case kJCluster:
            return fInputListCluster;
        case kJPtMatchedPi0:
            return fInputListPtMatchedPi0;
		default:
			std::cout << "Particle species not specified, return nothing" << std::endl;
			return 0;
	}
}

Int_t AliJHMRGeantCatalyst::LoadInput()
{
    fIn = TFile::Open(inputfile.Data(), "READ");
    if (!fIn) return 0;

    fTree = (TTree*)fIn->Get("Run");
    fTree->SetBranchAddress("Events", &fEvent);
    fBrtrack = fTree->GetBranch("Events.fTracks");
    fBrtrack->SetAddress(&fTracks);
    fBrcluster = fTree->GetBranch("Events.fClusters");
    fBrcluster->SetAddress(&fClusters);
    return 1;
}