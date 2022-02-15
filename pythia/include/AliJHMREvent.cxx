#include "AliJHMREvent.h"

AliJHMREvent::AliJHMREvent()
{
    fNtrack = 0;
    fNcluster = 0;
    fTracks = new TClonesArray("AliJHMRTrack", 1000);
    fClusters = new TClonesArray("AliJHMRCluster", 1000);
}

AliJHMREvent::~AliJHMREvent()
{
    Clear();
}

void AliJHMREvent::Clear(Option_t * /*option*/)
{
    fNtrack = 0;
    fNcluster = 0;
    fTracks->Clear("C");
    fClusters->Clear("C");
}

AliJHMRTrack * AliJHMREvent::AddTrack(Int_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t vx, Float_t vy, Float_t vz, Float_t charge, Float_t eta)
{
    AliJHMRTrack *track = (AliJHMRTrack*)fTracks->ConstructedAt(fNtrack++);
    track->Set(pid, px, py, pz, e, vx, vy, vz, charge, eta);
    return track;
}

AliJHMRCluster * AliJHMREvent::AddCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal)
{
    AliJHMRCluster *clust = (AliJHMRCluster*)fClusters->ConstructedAt(fNcluster++);
    clust->Set(x, y, z, e, ehcal);
    return clust;
}

ClassImp(AliJHMREvent)
