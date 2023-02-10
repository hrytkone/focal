#ifndef ALIJHMREVENT_H
#define ALIJHMREVENT_H

#include <iostream>

#include "TObject.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "AliJBaseTrack.h"
#include "AliJHMRCluster.h"

class AliJBaseTrack;
class AliJHMRCluster;

// -----------------
// |  Event class  |
// -----------------
class AliJHMREvent : public TObject {

    private:
        static const Int_t NSEGMENT = 10;

        Int_t fNtrack;
        Int_t fNcluster;
        Int_t UniqueID;
        TClonesArray *fTracks;
        TClonesArray *fClusters;

        AliJHMREvent(const AliJHMREvent&) = delete;
        AliJHMREvent &operator=(const AliJHMREvent&) = delete;

    public:
        AliJHMREvent();
        ~AliJHMREvent();
        void Clear(Option_t * option="");

        void InitEvent();
        void SetNtrack(Int_t n) { fNtrack = n; }
        void SetNcluster(Int_t n) { fNcluster = n; }
        Int_t GetNtrack() { return fNtrack; }
        Int_t GetNcluster() { return fNcluster; }
        AliJBaseTrack * AddTrack(Int_t pid, Float_t px, Float_t py, Float_t pz, Float_t e, Float_t charge, Int_t partIdx, Int_t motherId, Int_t motherType, bool isPrimary);
        AliJHMRCluster * AddCluster(Float_t x, Float_t y, Float_t z, Float_t e, Float_t ehcal, Float_t eseg[], Float_t eseed[], Float_t width1[], Float_t width2[]);

        ClassDef(AliJHMREvent, 2)
};

#endif
