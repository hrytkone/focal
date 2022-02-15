#include "AliJHMRTrack.h"

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
}

void AliJHMRTrack::Clear(Option_t * /*option*/)
{
    TObject::Clear();
}

ClassImp(AliJHMRTrack)
