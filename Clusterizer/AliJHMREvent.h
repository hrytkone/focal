#ifndef ALIJHMREVENT_H
#define ALIJHMREVENT_H

#include <iostream>

#include "TObject.h"
#include "TClonesArray.h"
#include "TMath.h"

// -----------------
// |  Track class  |
// -----------------

class AliJHMRTrack : public TObject {

    private:
        Float_t fPx;       // x component of the momentum
        Float_t fPy;       // y component of the momentum
        Float_t fPz;       // z component of the momentum
        Float_t fVx;       // x vertex position
        Float_t fVy;       // y vertex position
        Float_t fVz;       // z vertex position
        Float_t fCharge;   // charge of the track
        Float_t fEnergy;   // energy of the track
        Float_t fEta;      // pseudorapidity of the track
        Int_t fPid;        // particle id

    public:
        AliJHMRTrack() { fPx = 0; fPy = 0; fPz = 0; fVx = 0; fVy = 0; fVz = 0; fCharge = 0; fEnergy = 0; fPid = 0; fEta = 0; }
        AliJHMRTrack(Float_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t vx, Float_t vy, Float_t vz, Float_t charge, Float_t eta);
        virtual ~AliJHMRTrack() { Clear(); }
        void Clear(Option_t * option="");

        void Set(Float_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t vx, Float_t vy, Float_t vz, Float_t charge, Float_t eta);

        Float_t GetPx() const { return fPx; }
        Float_t GetPy() const { return fPy; }
        Float_t GetPz() const { return fPz; }
        Float_t GetVx() const { return fVx; }
        Float_t GetVy() const { return fVy; }
        Float_t GetVz() const { return fVz; }
        Float_t GetPt() const { return TMath::Sqrt(fPx*fPx + fPy*fPy); }
        Float_t GetCharge() const { return fCharge; }
        Float_t GetE() const { return fEnergy; }
        Float_t GetEta() const { return fEta; }
        Int_t GetPid() const { return fPid; }

        ClassDef(AliJHMRTrack, 1)
};

// -------------------
// |  Cluster class  |
// -------------------

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

// -----------------
// |  Event class  |
// -----------------
class AliJHMREvent : public TObject {

    private:
        Int_t fNtrack;
        Int_t fNcluster;
        TClonesArray *fTracks;
        TClonesArray *fClusters;

        AliJHMREvent(const AliJHMREvent&) = delete;
        AliJHMREvent &operator=(const AliJHMREvent&) = delete;

    public:
        AliJHMREvent();
        ~AliJHMREvent();
        void Clear(Option_t * option="");

        void SetNtrack(Int_t n) { fNtrack = n; }
        void SetNcluster(Int_t n) { fNcluster = n; }
        Int_t GetNtrack() { return fNtrack; }
        Int_t GetNcluster() { return fNcluster; }
        AliJHMRTrack * AddTrack(Int_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t vx, Float_t vy, Float_t vz, Float_t charge, Float_t eta);
        AliJHMRCluster * AddCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal);

        ClassDef(AliJHMREvent, 1)
};

#endif
