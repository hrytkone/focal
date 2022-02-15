#ifndef ALIJHMRTRACK_H
#define ALIJHMRTRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#include <TString.h>
#include <TMath.h>
#endif

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

#endif
