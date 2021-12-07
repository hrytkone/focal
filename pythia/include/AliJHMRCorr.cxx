#include "AliJHMRCorr.h"

bool AliJHMRCorr::IsTrackerAcceptance(double eta, double etaRange)
{
    return TMath::Abs(eta) < etaRange/2. ? true : false;
}

bool AliJHMRCorr::IsDetAcceptance(double eta, detector labelDet)
{
    return (eta > detEta[labelDet][0] && eta < detEta[labelDet][1]) ? true : false;
}

bool AliJHMRCorr::IsMassWindow(double mass)
{
    return (mass > 110. && mass < 160.) ? true : false;
}

bool AliJHMRCorr::IsSideband(double mass)
{
    return ((mass > 40. && mass < 80.) || (mass > 210. && mass < 280.)) ? true : false;
}

// Return 1 if trigger from peak has larger pT than trigger from sideband, 0 if sideband
// trigger pT is larger
int AliJHMRCorr::GetLargerTrigg(TClonesArray *arrPi0Peak, std::vector<int> listTriggPeak, TClonesArray *arrPi0Side, std::vector<int> listTriggSide)
{
    if (!listTriggPeak.size()) return 0;
    if (!listTriggSide.size()) return 1;
    TLorentzVector *lvTriggPeak = (TLorentzVector*)arrPi0Peak->At(listTriggPeak[0]);
    TLorentzVector *lvTriggSide = (TLorentzVector*)arrPi0Side->At(listTriggSide[0]);
    double pTpeak = lvTriggPeak->Pt();
    double pTside = lvTriggSide->Pt();
    if (pTpeak > pTside) return 1;
    return 0;
}

void AliJHMRCorr::ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, bool bMass)
{
    int nCandidate = 0;
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector lvSum = GetPhotonSumVector(arrPhoton, i, j);
            double mass = 1000.*lvSum.M();
            //double pT = lvSum.Pt();
            bool bIsInWindow = bMass ? IsMassWindow(mass) : IsSideband(mass);
            if (bIsInWindow) new ((*arrPi0Candidates)[nCandidate++]) TLorentzVector(lvSum);
        }
    }
}

void AliJHMRCorr::GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bUseLeading)
{
    int iLeadingTrigg = -1;
    if (bUseLeading) {
        iLeadingTrigg = GetLeadingTriggerIndex(arrPi0Candidates);
        if (iLeadingTrigg > -1) listTrigg.push_back(iLeadingTrigg);
    }

    int nPi0Candidates = arrPi0Candidates->GetEntriesFast();
    for (int i = 0; i < nPi0Candidates; i++) {
        TLorentzVector *lvPion = (TLorentzVector*)arrPi0Candidates->At(i);
        double pT = lvPion->Pt();
        int iTrigg = GetBin(triggPt, NTRIGGBINS, pT);
        int iAssoc = GetBin(assocPt, NASSOCBINS, pT);
        if (bUseLeading) {
            if (i != iLeadingTrigg && iAssoc >= 0) listAssoc.push_back(i);
        } else {
            if (iTrigg >= 0) listTrigg.push_back(i);
            if (iAssoc >= 0) listAssoc.push_back(i);
        }
        if (iTrigg > 0) binsWithTrigg[iTrigg]++;
    }
}

void AliJHMRCorr::DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bUseWeight)
{
    double wTrigg = 1.0;
    double wAssoc = 1.0;

    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrPi0->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();

        int iTriggBin = GetBin(triggPt, NTRIGGBINS, ptTrigg);

        int nAssoc = listAssoc.size();
        if (nAssoc < 1) continue;

        if (bUseWeight) wTrigg = 1./fPhotonAcceptanceEfficiency->Eval(ptTrigg);

        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];
            if (iTrigg==iAssoc) continue; // autocorrelations
            TLorentzVector *lvAssoc = (TLorentzVector*)arrPi0->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, NASSOCBINS, ptAssoc);

            if (triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;

            if (bUseWeight) wAssoc = 1./fPhotonAcceptanceEfficiency->Eval(ptAssoc);

            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;
            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
        }
    }

}

void AliJHMRCorr::DoCorrelations(TClonesArray *arrPi0Trigg, std::vector<int> listTrigg, TClonesArray *arrPi0Assoc, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bUseWeightTrigg, bool bUseWeightAssoc)
{
    double wTrigg = 1.0;
    double wAssoc = 1.0;

    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrPi0Trigg->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();
        int iTriggBin = GetBin(triggPt, NTRIGGBINS, ptTrigg);

        int nAssoc = listAssoc.size();
        if (nAssoc < 1) continue;

        if (bUseWeightTrigg) wTrigg = 1./fPhotonAcceptanceEfficiency->Eval(ptTrigg);

        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];

            TLorentzVector *lvAssoc = (TLorentzVector*)arrPi0Assoc->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, NASSOCBINS, ptAssoc);

            if (triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;

            if (bUseWeightAssoc) wAssoc = 1./fPhotonAcceptanceEfficiency->Eval(ptAssoc);

            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;
            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
        }
    }
}

int AliJHMRCorr::GetBin(double arr[], int nArr, double val)
{
    for (int i=0; i<nArr; i++) {
        if (arr[i]<=val && val<arr[i+1]) return i;
    }
    return -1;
}

double AliJHMRCorr::GetDeltaPhi(double phiTrigg, double phiAssoc)
{
    double dphi = phiTrigg - phiAssoc;
    if (dphi>3.0*TMath::Pi()/2.0) return dphi - 2.0*TMath::Pi();
    if (dphi<-TMath::Pi()/2.0) return dphi + 2.0*TMath::Pi();
    return dphi;
}

// For FoCal:
//      sigma/E = a/sqrt(E) + b = 0.27/sqrt(E) + 0.01
double AliJHMRCorr::PhotonEnergySmearing(double px, double py, double pz)
{
    const double eCut = 0.01;
    double e = sqrt(px*px + py*py + pz*pz);
    if (e < eCut) return e;

    double a = 0.27;
    double b = 0.01;
    double sigma = sqrt(a*a*e + b*b*e*e);
    double eSmear = -1.0;
    while (eSmear < 0.0) eSmear = fRand->Gaus(e, sigma);
    return eSmear;
}

void AliJHMRCorr::SmearEnergies(TClonesArray * arrParticles)
{
    int nent = arrParticles->GetEntriesFast();
    for ( int ient = 0; ient < nent; ient++ ) {
        TLorentzVector *lv = (TLorentzVector*)arrParticles->At(ient);
        double px = lv->Px();
        double py = lv->Py();
        double pz = lv->Pz();
        double e = lv->E();
        double eSmear = PhotonEnergySmearing(px, py, pz);
        lv->SetPxPyPzE(eSmear*px/e, eSmear*py/e, eSmear*pz/e, eSmear);
    }
}

// To count in single photon efficiency
bool AliJHMRCorr::IsPhotonRemoved(double ePhoton)
{
    double photonEff = fPhotonEfficiency->Eval(ePhoton);
    if (fRand->Uniform() > photonEff) return true;
    return false;
}

TLorentzVector AliJHMRCorr::GetPhotonSumVector(TClonesArray *arrPhoton, int iPhoton1, int iPhoton2)
{
    TLorentzVector *lv1 = (TLorentzVector*)arrPhoton->At(iPhoton1);
    TLorentzVector *lv2 = (TLorentzVector*)arrPhoton->At(iPhoton2);
    TLorentzVector lvSum;
    lvSum = ((*lv1) + (*lv2));
    return lvSum;
}

int AliJHMRCorr::GetLeadingTriggerIndex(TClonesArray *arrPi0)
{
    int itrigg = -1;
    double pTmax = 0;
    int nPi0 = arrPi0->GetEntriesFast();
    for (int i = 0; i < nPi0; i++) {
        TLorentzVector *lvPi0 = (TLorentzVector*)arrPi0->At(i);
        double pT = lvPi0->Pt();
        if (pT > triggPt[0] && pTmax < pT) {
            pTmax = pT;
            itrigg = i;
        }
    }
    return itrigg;
}

void AliJHMRCorr::FillRealTriggers(TClonesArray *arrRealPi0, std::vector<int>& listTrigg)
{
    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        TLorentzVector *lvTrigg = (TLorentzVector*)arrRealPi0->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        int iTriggBin = GetBin(triggPt, NTRIGGBINS, ptTrigg);
        histos->hRealTriggCounter->Fill(iTriggBin+0.5);
    }
}

void AliJHMRCorr::FillPionMasses(TClonesArray *arrPhoton, int binsWithTriggPeak[NTRIGGBINS], int binsWithTriggSide[NTRIGGBINS])
{
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        for (int j = 0; j < i; j++) {
            TLorentzVector lvSum = GetPhotonSumVector(arrPhoton, i, j);
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            int iTriggBin = GetBin(triggPt, NTRIGGBINS, pT);
            int iAssocBin = GetBin(assocPt, NASSOCBINS, pT);
            if (iTriggBin >= 0) histos->hPi0MassTrigg[iTriggBin]->Fill(mass);
            for (int it = 0; it < NTRIGGBINS; it++) {
                if (triggPt[it] < assocPt[iAssocBin+1]) continue;
                if (binsWithTriggPeak[it] > 0 && iAssocBin >= 0) histos->hPi0MassAssocPeak[it][iAssocBin]->Fill(mass);
                if (binsWithTriggSide[it] > 0 && iAssocBin >= 0) histos->hPi0MassAssocSide[it][iAssocBin]->Fill(mass);
            }
        }
    }
}
