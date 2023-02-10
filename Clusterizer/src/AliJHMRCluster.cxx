#include "AliJHMRCluster.h"

#ifdef __MAKECINT__
#pragma link C++ class AliJHMRCluster;
#endif

AliJHMRCluster::AliJHMRCluster() {
    fX = 0; fY = 0; fZ = 0; fEnergy = 0; fEnergyHCAL = 0;
    for (Int_t i = 0; i < NSEGMENT; i++) {
        fSegmentEnergy[i] = 0;
        fSeedEnergy[i] = 0;
        fWidth1[i] = 0;
        fWidth2[i] = 0;
    }
}

AliJHMRCluster::AliJHMRCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal)
{
    fX = x;
    fY = y;
    fZ = z;
    fEnergy = e;
    fEnergyHCAL = ehcal;
    for (Int_t i = 0; i < NSEGMENT; i++) {
        fSegmentEnergy[i] = 0;
        fSeedEnergy[i] = 0;
        fWidth1[i] = 0;
        fWidth2[i] = 0;
    }
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

void AliJHMRCluster::SetSegmentEnergy(Int_t iseg, Float_t e)
{
    fSegmentEnergy[iseg] = e;
}

void AliJHMRCluster::SetSeedEnergy(Int_t iseg, Float_t e)
{
    fSeedEnergy[iseg] = e;
}

void AliJHMRCluster::SetWidth1(Int_t iseg, Float_t width1)
{
    fWidth1[iseg] = width1;
}

void AliJHMRCluster::SetWidth2(Int_t iseg, Float_t width2)
{
    fWidth2[iseg] = width2;
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
