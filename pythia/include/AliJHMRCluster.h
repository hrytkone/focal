#ifndef ALIJHMRCLUSTER_H
#define ALIJHMRCLUSTER_H

#ifndef ROOT_TObject
#include <TObject.h>
#include <TString.h>
#include <TMath.h>
#endif

class AliJHMRCluster : public TObject {

    private:
        Float_t fX;            // cluster x component
        Float_t fY;            // cluster y component
        Float_t fZ;            // cluter z component
        Float_t fEnergy;       // cluster energy in ECAL
        Float_t fEnergyHCAL;   // cluster energy in HCAL

    public:
        AliJHMRCluster() { fX = 0; fY = 0; fZ = 0; fEnergy = 0; fEnergyHCAL = 0; }
        AliJHMRCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal);

        virtual ~AliJHMRCluster() { Clear(); }

        void Clear(Option_t * option="");

        void Set(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal);

        Float_t GetX() const { return fX; }
        Float_t GetY() const { return fY; }
        Float_t GetZ() const { return fZ; }
        Float_t GetE() const { return fEnergy; }
        Float_t GetEHCAL() const { return fEnergyHCAL; }
        Float_t GetPhi() const { return TMath::ATan2(fY, fX); }
        Float_t GetPt() const;
        Float_t GetEta() const;

        ClassDef(AliJHMRCluster, 1)
};

#endif
