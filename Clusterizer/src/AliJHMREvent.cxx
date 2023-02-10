#include "AliJHMREvent.h"

#ifdef __MAKECINT__
//#pragma link C++ class AliJBaseTrack;
//#pragma link C++ class AliJHMRCluster;
#pragma link C++ class AliJHMREvent+;
#endif

//_____________________________________________________________________________________
// Event class methods
//

AliJHMREvent::AliJHMREvent()
{
    fNtrack = 0;
    fNcluster = 0;
    UniqueID = 0;
    fTracks = 0;
    fClusters = 0;
}

AliJHMREvent::~AliJHMREvent()
{
    Clear();
}

void AliJHMREvent::InitEvent()
{
    fTracks = new TClonesArray("AliJBaseTrack", 1000);
    fClusters = new TClonesArray("AliJHMRCluster", 1000);
}

void AliJHMREvent::Clear(Option_t * /*option*/)
{
    fNtrack = 0;
    fNcluster = 0;
    fTracks->Clear("C");
    fClusters->Clear("C");
}

AliJBaseTrack *AliJHMREvent::AddTrack(Int_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t charge, Int_t partIdx, Int_t motherId, Int_t motherType, bool isPrimary)
{
    AliJBaseTrack *track = (AliJBaseTrack*)fTracks->ConstructedAt(fNtrack++);
    track->SetPxPyPzE(px, py, pz, e);
    track->SetUniqueID(UniqueID++);
    track->SetID(pid);
    track->SetMCIndex(partIdx);
    track->SetMotherID(motherId);
    track->SetMotherType(motherType);
    track->SetCharge(charge);
    track->SetPrimary(isPrimary);
    track->SetTrackEff(1.);
    return track;
}

AliJHMRCluster * AliJHMREvent::AddCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal, Float_t eseg[], Float_t eseed[], Float_t width1[], Float_t width2[])
{
    AliJHMRCluster *clust = (AliJHMRCluster*)fClusters->ConstructedAt(fNcluster++);
    clust->Set(x, y, z, e, ehcal);
    for (Int_t i=0; i<NSEGMENT; i++) {
        clust->SetSegmentEnergy(i, eseg[i]);
        clust->SetSeedEnergy(i, eseed[i]);
        clust->SetWidth1(i, width1[i]);
        clust->SetWidth2(i, width2[i]);
    }
    return clust;
}
