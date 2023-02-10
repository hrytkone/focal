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

void AliJHMRGeantCatalyst::GetPtMatchedClusters()
{
    std::vector<bool> bParticleIsMatched;
    for (int i = 0; i < fInputListPhoton->GetEntriesFast(); i++) bParticleIsMatched.push_back(false);
    //cout << "NCLUSTER : " << fInputListCluster->GetEntriesFast() << "  NPHOTON : " << fInputListPhoton->GetEntriesFast() << endl;
    for (int icluster = 0; icluster < fInputListCluster->GetEntriesFast(); icluster++) {
        AliJBaseTrack *cluster = (AliJBaseTrack*)fInputListCluster->At(icluster);
        double etaCluster = cluster->Eta();
        double phiCluster = cluster->Phi();

        // Get geometrically closest MC photon to the cluster
        double dist = 1000000.0;
        double ptMatched = -1;
        int iMatched = 0;
        for (int itr = 0; itr < fInputListPhoton->GetEntriesFast(); itr++) {
            AliJBaseTrack *tr = (AliJBaseTrack*)fInputListPhoton->At(itr);

            double etaPhoton = tr->Eta();
            double phiPhoton = tr->Phi();
            double ptPhoton = tr->Pt();
            double deltaPhoton = TMath::Sqrt((etaPhoton-etaCluster)*(etaPhoton-etaCluster) + (phiPhoton-phiCluster)*(phiPhoton-phiCluster));

            if (deltaPhoton < dist) {
                dist = deltaPhoton;
                ptMatched = ptPhoton;
                if (!bParticleIsMatched[iMatched]) iMatched = itr;
            }
        }
        bParticleIsMatched[iMatched] = true;

        AliJBaseTrack *matched = (AliJBaseTrack*)fInputListPhoton->At(iMatched);
        //cout << "Cluster (phi eta pt): (" << phiCluster << " " << etaCluster << " " << cluster->Pt() << ")\t particle (phi eta pt): (" << matched->Phi() << " " << matched->Eta() << " " << ptMatched << ")\t dist = " << dist << endl;

        // Save cluster to the list but with matched pt
        lvParticle.SetPtEtaPhiM(ptMatched, cluster->Eta(), cluster->Phi(), 0.);
        lvParticle.SetUniqueID(UniqueID++);

        AliJBaseTrack track( lvParticle );
        new((*fInputListPtMatchedCluster)[fInputListPtMatchedCluster->GetEntriesFast()]) AliJBaseTrack(track);        
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
        case kJPtMatchedCluster:
            return fInputListPtMatchedCluster;
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