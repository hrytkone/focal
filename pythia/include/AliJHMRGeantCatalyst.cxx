#include "AliJHMRGeantCatalyst.h"
#include <set>

using namespace std;

TClonesArray * AliJHMRGeantCatalyst::GetClusters()
{
    for (int iclust = 0; iclust < fClusters->GetEntries(); iclust++) {
        AliJHMRCluster *cluster = (AliJHMRCluster*)fClusters->At(iclust);
        lvParticle.SetPtEtaPhiE(cluster->GetPt(), cluster->GetEta(), cluster->GetPhi(), cluster->GetE());
        lvParticle.SetUniqueID(UniqueID++);

        AliJBaseTrack track( lvParticle );
        new((*fInputListCluster)[fInputListCluster->GetEntriesFast()]) AliJBaseTrack(track);
    }
    return fInputListCluster;
}

TClonesArray * AliJHMRGeantCatalyst::GetPi0True()
{
    for (int itrack = 0; itrack < fTracks->GetEntries(); itrack++) {
        AliJHMRTrack *tr = (AliJHMRTrack*)fTracks->At(itrack);
        lvParticle.SetPxPyPzE(tr->GetPx(), tr->GetPy(), tr->GetPz(), tr->GetE());
        lvParticle.SetUniqueID(UniqueID++);

        AliJBaseTrack track( lvParticle );
        new((*fInputListTrack)[fInputListTrack->GetEntriesFast()]) AliJBaseTrack(track);
    }
    return fInputListTrack;
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
