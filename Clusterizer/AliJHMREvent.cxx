#include "AliJHMREvent.h"

#ifdef __MAKECINT__
#pragma link C++ class AliJHMREvent+;
#pragma link C++ class AliJHMRTrack+;
#pragma link C++ class AliJHMRCluster+;
#endif

//_____________________________________________________________________________________
// Event class methods
//

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

//_____________________________________________________________________________________
// Track class methods
//

AliJHMRTrack::AliJHMRTrack(Float_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t vx, Float_t vy, Float_t vz, Float_t charge, Float_t eta)
{
    fPid = pid;
    fPx = px;
    fPy = py;
    fPz = pz;
    fEnergy = e;
    fVx = vx;
    fVy = vy;
    fVz = vz;
    fCharge = charge;
    fEta = eta;
}

void AliJHMRTrack::Set(Float_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t vx, Float_t vy, Float_t vz, Float_t charge, Float_t eta)
{
    fPid = pid;
    fPx = px;
    fPy = py;
    fPz = pz;
    fEnergy = e;
    fVx = vx;
    fVy = vy;
    fVz = vz;
    fCharge = charge;
    fEta = eta;
}

void AliJHMRTrack::Clear(Option_t * /*option*/)
{
    TObject::Clear();
}

//_____________________________________________________________________________________
// Cluster class methods
//

AliJHMRCluster::AliJHMRCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal)
{
    fX = x;
    fY = y;
    fZ = z;
    fEnergy = e;
    fEnergyHCAL = ehcal;
}

void AliJHMRCluster::Clear(Option_t * /*option*/)
{
    TObject::Clear();
}

void AliJHMRCluster::Set(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal)
{
    fX = x;
    fY = y;
    fZ = z;
    fEnergy = e;
    fEnergyHCAL = ehcal;
}

Float_t AliJHMRCluster::GetPt() const
{
    return (fEnergy/TMath::Sqrt(fX*fX + fY*fY + fZ*fZ)) * TMath::Sqrt(fX*fX + fY*fY);
}

Float_t AliJHMRCluster::GetEta() const
{
    float theta = 0.;
    if (fZ != 0) theta = TMath::ATan(TMath::Sqrt(fX*fX + fY*fY)/fZ);
    if (theta > 1e-6)
        return -TMath::Log(TMath::Tan(theta/2.));
    else
        return -1000.;
}
