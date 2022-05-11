#ifndef ALIJHMRCLUSTER_H
#define ALIJHMRCLUSTER_H

#include <TObject.h>
#include <TString.h>
#include <TMath.h>

class AliJHMRCluster : public TObject {

    private:

        static const Int_t NSEGMENT = 10;

        Float_t fX;            // cluster x component
        Float_t fY;            // cluster y component
        Float_t fZ;            // cluter z component
        Float_t fEnergy;       // cluster energy in ECAL
        Float_t fEnergyHCAL;   // cluster energy in HCAL

        Float_t fSegmentEnergy[NSEGMENT]; // cluster energy per segment
        Float_t fSeedEnergy[NSEGMENT];    // seed energy per segment
        Float_t fWidth1[NSEGMENT];        // major axis width of the cluster
        Float_t fWidth2[NSEGMENT];        // minor axis width of the cluster

    public:
        AliJHMRCluster();
        AliJHMRCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal);

        virtual ~AliJHMRCluster() { Clear(); }

        void Clear(Option_t * option="");

        void Set(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal);
        void SetSegmentEnergy(Int_t iseg, Float_t e);
        void SetSeedEnergy(Int_t iseg, Float_t e);
        void SetWidth1(Int_t iseg, Float_t width1);
        void SetWidth2(Int_t iseg, Float_t width2);

        Float_t GetSegmentEnergy(Int_t iseg) const { return fSegmentEnergy[iseg]; }
        Float_t GetSeedEnergy(Int_t iseg) const { return fSeedEnergy[iseg]; }
        Float_t GetWidth1(Int_t iseg) const { return fWidth1[iseg]; }
        Float_t GetWidth2(Int_t iseg) const { return fWidth2[iseg]; }

        Float_t GetX() const { return fX; }
        Float_t GetY() const { return fY; }
        Float_t GetZ() const { return fZ; }
        Float_t GetE() const { return fEnergy; }
        Float_t GetEHCAL() const { return fEnergyHCAL; }
        Float_t GetPhi() const { return TMath::ATan2(fY, fX); }
        Float_t GetPt() const;
        Float_t GetEta() const;

        ClassDef(AliJHMRCluster, 2)
};

#endif
