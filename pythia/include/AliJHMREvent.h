#ifndef ALIJHMREVENT_H
#define ALIJHMREVENT_H

#include <iostream>

#include "TObject.h"
#include "TClonesArray.h"
#include "TMath.h"

#include "AliJHMRCluster.h"
#include "AliJHMRTrack.h"

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
