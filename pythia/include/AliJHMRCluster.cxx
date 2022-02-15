#include "AliJHMRCluster.h"

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

ClassImp(AliJHMRCluster)
