#include "AliJHMRCorr.h"

bool AliJHMRCorr::IsTrackerAcceptance(double eta, double etaRange)
{
    return TMath::Abs(eta) < etaRange/2. ? true : false;
}

bool AliJHMRCorr::IsDetAcceptance(double eta, detector labelDet)
{
    return (eta > detEta[labelDet][0] && eta < detEta[labelDet][1]) ? true : false;
}

// Fixed mass window (1. check for full sim)
// done in reconstruction phase
bool AliJHMRCorr::IsMassWindow(double mass)
{
    return (mass > massWindowMin && mass < massWindowMax) ? true : false;
}

// 2. check for full sim, use mass window of 3sigmas (in full sim), different for each pt bin
// done in trigger-assoc division phase
bool AliJHMRCorr::IsMassWindow(double mass, int ibin, bool isTriggBin)
{
    if (isTriggBin) {
        if (mass > massPeakPosTrigg[ibin]-3.*massSigmaTrigg[ibin] && mass < massPeakPosTrigg[ibin]+3.*massSigmaTrigg[ibin])
            return true;
        else
            return false;
    } else {
        if (mass > massPeakPosAssoc[ibin]-3.*massSigmaAssoc[ibin] && mass < massPeakPosAssoc[ibin]+3.*massSigmaAssoc[ibin])
            return true;
        else
            return false;
    }
    return false;
}

bool AliJHMRCorr::IsSideband(double mass)
{
    //if (fIsFullSim)
        return (mass > sidebandMin && mass < sidebandMax) ? true : false;
    //else
    //    return ((mass > 40. && mass < 80.) || (mass > 210. && mass < 280.)) ? true : false;
}

// Return 1 if trigger from peak has larger pT than trigger from sideband, 0 if sideband
// trigger pT is larger
int AliJHMRCorr::GetLargerTrigg(TClonesArray *arrPi0Peak, std::vector<int> listTriggPeak, TClonesArray *arrPi0Side, std::vector<int> listTriggSide)
{
    if (!listTriggPeak.size() && !listTriggSide.size()) return -1; // no trigger in either
    if (!listTriggPeak.size()) return 0;
    if (!listTriggSide.size()) return 1;
    AliJBaseTrack *lvTriggPeak = (AliJBaseTrack*)arrPi0Peak->At(listTriggPeak[0]);
    AliJBaseTrack *lvTriggSide = (AliJBaseTrack*)arrPi0Side->At(listTriggSide[0]);
    double pTpeak = lvTriggPeak->Pt();
    double pTside = lvTriggSide->Pt();
    if (pTpeak > pTside) return 1;
    return 0;
}

int AliJHMRCorr::ReconstructPions(TClonesArray *arrPhoton, TClonesArray *arrPi0Candidates, detector idet, bool bMass)
{
    int nTrue = 0;
    int nCandidate = 0;
    int nPhoton = arrPhoton->GetEntriesFast();

    for (int i = 1; i < nPhoton; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)arrPhoton->At(i);
        if (!fIsFullSim && IsPhotonRemoved(lv1->E())) continue;
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)arrPhoton->At(j);
            if (!fIsFullSim && IsPhotonRemoved(lv2->E())) continue;
            if (GetAsymmetry(arrPhoton, lv1, lv2)>asymcut) continue;
            AliJBaseTrack lvSum = GetPhotonSumVector(arrPhoton, lv1, lv2);
            if (lvSum.Eta()<detEta[idet][0]+etacut || lvSum.Eta()>detEta[idet][1]-etacut) continue;
            double mass = 1000.*lvSum.M();
            bool bIsInWindow = bMass ? IsMassWindow(mass) : IsSideband(mass);
            if (bIsInWindow) {
                new ((*arrPi0Candidates)[nCandidate++]) AliJBaseTrack(lvSum);
                if (lvSum.GetLabel()==1) nTrue++;
                std::vector<int> pair = {i,j};
                if (bMass)
                    photonId.push_back(pair);
                else
                    sidebandId.push_back(pair);
            }
        }
    }
    return nTrue;
}

void AliJHMRCorr::GetTriggAssocLists(TClonesArray *arrPi0Candidates, std::vector<int>& listTrigg, std::vector<int>& listAssoc, int *binsWithTrigg, bool bMass)
{
    int iLeadingTrigg = -1;
    if (fUseLeading) {
        iLeadingTrigg = GetLeadingTriggerIndex(arrPi0Candidates, bMass);
        if (iLeadingTrigg > -1) {
            listTrigg.push_back(iLeadingTrigg);
            binsWithTrigg[0]++;
        }
    }

    int nPi0Candidates = arrPi0Candidates->GetEntriesFast();
    for (int i = 0; i < nPi0Candidates; i++) {
        AliJBaseTrack *lvPion = (AliJBaseTrack*)arrPi0Candidates->At(i);
        double pT = lvPion->Pt();
        double mass = 1000.*lvPion->M();
        int iTrigg = GetBin(triggPt, NTRIGGBINS, pT);
        int iAssoc = GetBin(assocPt, NASSOCBINS, pT);
        if (fUseLeading) {
            if (fIsFullSim && bMass) {
                int ibin = GetBin(leadingPt, NLEADINGBINS, pT);
                if (i != iLeadingTrigg && iAssoc >= 0 && IsMassWindow(mass, ibin, 1)) listAssoc.push_back(i);
            } else {
                if (i != iLeadingTrigg && iAssoc >= 0) listAssoc.push_back(i);
            }
        } else {
            if (fIsFullSim && bMass) {
                if (iTrigg >= 0 && IsMassWindow(mass, iTrigg, 1)) listTrigg.push_back(i);
                if (iAssoc >= 0 && IsMassWindow(mass, iAssoc, 0)) listAssoc.push_back(i);
            } else {
                if (iTrigg >= 0) listTrigg.push_back(i);
                if (iAssoc >= 0) listAssoc.push_back(i);
            }
            if (bMass) {
                if (iTrigg >= 0 && IsMassWindow(mass, iTrigg, 1)) binsWithTrigg[iTrigg]++;
            } else {
                if (iTrigg >= 0) binsWithTrigg[iTrigg]++;
            }
        }
    }
}

void AliJHMRCorr::DoCorrelations(TClonesArray *arrPi0, std::vector<int> listTrigg, std::vector<int> listAssoc, TH2D *hCorr[NTRIGGBINS][NASSOCBINS], bool bTrueCorr, bool bMassWindowTrigg, bool bUseWeight)
{
    double wTrigg = 1.0;
    double wAssoc = 1.0;

    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    bool etriggFilled = false;
    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        AliJBaseTrack *lvTrigg = (AliJBaseTrack*)arrPi0->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();

        int iTriggBin = GetBin(triggPt, NTRIGGBINS, ptTrigg);

        int nAssoc = listAssoc.size();
        if (nAssoc < 1) continue;

        //if (bUseWeight) wTrigg = 1./fPhotonAcceptanceEfficiency->Eval(ptTrigg);
        //if (bUseWeight && !fIsFullSim) wTrigg = 1./pi0eff;
        //if (bUseWeight && fIsFullSim) wTrigg = 1./effCorrTrigg[iTriggBin];

        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];
            if (iTrigg==iAssoc) continue; // autocorrelations
            AliJBaseTrack *lvAssoc = (AliJBaseTrack*)arrPi0->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, NASSOCBINS, ptAssoc);
            int iLeadingBin = GetBin(leadingPt, NLEADINGBINS, ptAssoc);

            if (!fUseLeading && triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;

            //if (bUseWeight) wAssoc = 1./fPhotonAcceptanceEfficiency->Eval(ptAssoc);
            if (bUseWeight && !fIsFullSim) wAssoc = 1./pi0eff;
            if (bUseWeight && fIsFullSim) {
                if (fUseLeading)
                    wAssoc = 1./effCorrLeading[iLeadingBin];
                else
                    wAssoc = 1./effCorrAssoc[iAssocBin];
            }

            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;

            if (CheckAssocPhotonPair(iTrigg, iAssoc, bMassWindowTrigg) || bTrueCorr)
                hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
            else if (bMassWindowTrigg)
                histos->hCorrRejectMassMass[iTriggBin][iAssocBin]->Fill(dphi, deta);
            else
                histos->hCorrRejectSideSide[iTriggBin][iAssocBin]->Fill(dphi, deta);

            // Histogram for checking the energies of sideband candidates
            if (!bTrueCorr && !bMassWindowTrigg) {
                if (!etriggFilled) {
                    histos->hEnergySidebandTrigg[iTriggBin]->Fill(lvTrigg->E());
                    etriggFilled = true;
                }
                histos->hEnergySidebandAssoc[iTriggBin][iAssocBin]->Fill(lvAssoc->E());
            }
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
        AliJBaseTrack *lvTrigg = (AliJBaseTrack*)arrPi0Trigg->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();
        int iTriggBin = GetBin(triggPt, NTRIGGBINS, ptTrigg);

        int nAssoc = listAssoc.size();
        if (nAssoc < 1) continue;

        //if (bUseWeightTrigg) wTrigg = 1./fPhotonAcceptanceEfficiency->Eval(ptTrigg);
        if (bUseWeightTrigg && !fIsFullSim) wTrigg = 1./pi0eff;
        if (bUseWeightTrigg && fIsFullSim) wTrigg = 1./effCorrTrigg[iTriggBin];

        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];

            AliJBaseTrack *lvAssoc = (AliJBaseTrack*)arrPi0Assoc->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, NASSOCBINS, ptAssoc);
            int iLeadingBin = GetBin(leadingPt, NLEADINGBINS, ptAssoc);

            if (!fUseLeading && triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;

            //if (bUseWeightAssoc) wAssoc = 1./fPhotonAcceptanceEfficiency->Eval(ptAssoc);
            if (bUseWeightAssoc && !fIsFullSim) wAssoc = 1./pi0eff;
            if (bUseWeightAssoc && fIsFullSim) {
                if (fUseLeading)
                    wAssoc = 1./effCorrLeading[iLeadingBin];
                else
                    wAssoc = 1./effCorrAssoc[iAssocBin];
            }

            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;
            hCorr[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
        }
    }
}

void AliJHMRCorr::ConstructTrueCorrComponents(TClonesArray *arrPi0, TClonesArray *arrPhoton, std::vector<int> listTrigg, std::vector<int> listAssoc, bool bUseWeight)
{
    double wTrigg = 1.0;
    double wAssoc = 1.0;

    int nTrigg = listTrigg.size();
    if (nTrigg < 1) return;

    for (int it = 0; it < nTrigg; it++) {
        int iTrigg = listTrigg[it];
        AliJBaseTrack *lvTrigg = (AliJBaseTrack*)arrPi0->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        double phiTrigg = lvTrigg->Phi();
        double etaTrigg = lvTrigg->Eta();

        int iTriggBin = GetBin(triggPt, NTRIGGBINS, ptTrigg);

        int nAssoc = listAssoc.size();
        if (nAssoc < 1) continue;

        //if (bUseWeight) wTrigg = 1./fPhotonAcceptanceEfficiency->Eval(ptTrigg);
        if (bUseWeight) wTrigg = 1./pi0eff;

        bool etriggFilled = false;
        for (int ia = 0; ia < nAssoc; ia++) {
            int iAssoc = listAssoc[ia];
            if (iTrigg==iAssoc) continue; // autocorrelations
            AliJBaseTrack *lvAssoc = (AliJBaseTrack*)arrPi0->At(iAssoc);
            double ptAssoc = lvAssoc->Pt();
            double phiAssoc = lvAssoc->Phi();
            double etaAssoc = lvAssoc->Eta();
            int iAssocBin = GetBin(assocPt, NASSOCBINS, ptAssoc);

            if (!fUseLeading && triggPt[iTriggBin] < assocPt[iAssocBin+1]) continue;

            //if (bUseWeight) wAssoc = 1./fPhotonAcceptanceEfficiency->Eval(ptAssoc);
            if (bUseWeight) wAssoc = 1./pi0eff;

            double dphi = GetDeltaPhi(phiTrigg, phiAssoc);
            double deta = etaTrigg - etaAssoc;

            if (!CheckAssocPhotonPair(iTrigg, iAssoc, 1)) continue;

            if (lvTrigg->GetLabel()==1 && lvAssoc->GetLabel()==1)
                histos->hCorrSignalSignal[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
            if (lvTrigg->GetLabel()==1 && lvAssoc->GetLabel()==0)
                histos->hCorrSignalBg[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
            if (lvTrigg->GetLabel()==0 && lvAssoc->GetLabel()==1)
                histos->hCorrBgSignal[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
            if (lvTrigg->GetLabel()==0 && lvAssoc->GetLabel()==0) {
                histos->hCorrBgBg[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
                if (!etriggFilled) {
                    histos->hEnergyMassBgTrigg[iTriggBin]->Fill(lvTrigg->E());
                    etriggFilled = true;
                }
                histos->hEnergyMassBgAssoc[iTriggBin][iAssocBin]->Fill(lvAssoc->E());

                std::vector<int> triggPairs = photonId[iTrigg];
                std::vector<int> assocPairs = photonId[iAssoc];
                AliJBaseTrack *t1 = (AliJBaseTrack*)arrPhoton->At(triggPairs[0]);
                AliJBaseTrack *t2 = (AliJBaseTrack*)arrPhoton->At(triggPairs[1]);
                AliJBaseTrack *a1 = (AliJBaseTrack*)arrPhoton->At(assocPairs[0]);
                AliJBaseTrack *a2 = (AliJBaseTrack*)arrPhoton->At(assocPairs[1]);

                if (a1->GetMotherID()==a2->GetMotherID()) {
                    cout << "Should not happen!" << endl;
                } else if (a1->GetLabel()==1 && a2->GetLabel()==1 && a1->GetMotherID()!=a2->GetMotherID()) {
                    histos->hCorrBgBgDecay[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
                    cout << "Mothers : " << endl;
                    cou t<< "\tt1=" << t1->GetMotherID() << "\tt2=" << t2->GetMotherID() << endl;
                    cou t<< "\ta1=" << a1->GetMotherID() << "\ta2=" << a2->GetMotherID() << endl;
                } else if (a1->GetLabel()!=a2->GetLabel()) {
                    histos->hCorrBgBgMix[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
                } else {
                    histos->hCorrBgBgNotDecay[iTriggBin][iAssocBin]->Fill(dphi, deta, wTrigg*wAssoc);
                }

            }
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

// since phi-acceptance is peridic
//double AliJHMRCorr::GetDeltaPhi(double phiTrigg, double phiAssoc)
//{
//    // dphi
//    double res = atan2(sin(phiTrigg-phiAssoc), cos(phiTrigg-phiAssoc));
//    return res > -1.413717 ? res : 6.283185+res ; // -9/20*pi : 2*pi
//}

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
        AliJBaseTrack *lv = (AliJBaseTrack*)arrParticles->At(ient);
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

// Cuts to energy & pT
bool AliJHMRCorr::IsPhotonRemoved(TClonesArray *arrPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    if (lv1->E() < 1. || lv2->E() < 1.) {
        return false;
    } else if (lv1->Pt() < 0.1 || lv2->Pt() < 0.1) {
        return false;
    } else {
        return true;
    }
}

double AliJHMRCorr::GetAsymmetry(TClonesArray *arrPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    return TMath::Abs(lv1->E() - lv2->E())/(lv1->E() + lv2->E());
}

AliJBaseTrack AliJHMRCorr::GetPhotonSumVector(TClonesArray *arrPhoton, AliJBaseTrack *lv1, AliJBaseTrack *lv2)
{
    AliJBaseTrack lvSum;
    lvSum = ((*lv1) + (*lv2));

    // Tag pi0 candidate with 1 if both photons in the pair are from same pi0 decay
    if (lv1->GetLabel()==1 && lv1->GetLabel()==1
        && lv1->GetMotherID()==lv2->GetMotherID())
        lvSum.SetLabel(1);
    else
        lvSum.SetLabel(0);

    return lvSum;
}

int AliJHMRCorr::GetLeadingTriggerIndex(TClonesArray *arrPi0, bool bUseSim)
{
    int itrigg = -1;
    double pTmax = 0;
    int nPi0 = arrPi0->GetEntriesFast();
    for (int i = 0; i < nPi0; i++) {
        AliJBaseTrack *lvPi0 = (AliJBaseTrack*)arrPi0->At(i);
        double pT = lvPi0->Pt();
        double mass = lvPi0->M();
        int iLeading = GetBin(leadingPt, NLEADINGBINS, pT);
        if (pT>leadingPt[NLEADINGBINS-1]) iLeading = NLEADINGBINS-1;
        if (bUseSim && !IsMassWindow(mass, iLeading, 1)) continue;
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
        AliJBaseTrack *lvTrigg = (AliJBaseTrack*)arrRealPi0->At(iTrigg);
        double ptTrigg = lvTrigg->Pt();
        int iTriggBin = GetBin(triggPt, NTRIGGBINS, ptTrigg);
        histos->hRealTriggCounter->Fill(iTriggBin+0.5);
    }
}

void AliJHMRCorr::FillPionMasses(TClonesArray *arrPhoton, int binsWithTriggPeak[NTRIGGBINS], int binsWithTriggSide[NTRIGGBINS], detector idet)
{
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)arrPhoton->At(i);
        if (!fIsFullSim && IsPhotonRemoved(lv1->E())) continue;
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)arrPhoton->At(j);
            if (GetAsymmetry(arrPhoton, lv1, lv2)>asymcut) continue;
            if (!fIsFullSim && IsPhotonRemoved(lv2->E())) continue;
            AliJBaseTrack lvSum = GetPhotonSumVector(arrPhoton, lv1, lv2);
            if (lvSum.Eta()<detEta[idet][0]+etacut || lvSum.Eta()>detEta[idet][1]-etacut) continue;
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();

            // Fill different components for the invariant mass distribution
            if (lv1->GetMotherID()==lv2->GetMotherID()) {
                histos->hMassTrue->Fill(mass);
            } else if (lv1->GetLabel()==1 && lv2->GetLabel()==1 && lv1->GetMotherID()!=lv2->GetMotherID()) {
                histos->hMassDecay->Fill(mass);
            } else if (lv1->GetLabel()!=lv2->GetLabel()) {
                histos->hMassMix->Fill(mass);
            } else {
                histos->hMassNotDecay->Fill(mass);
            }

            int iTriggBin = GetBin(triggPt, NTRIGGBINS, pT);
            int iAssocBin = GetBin(assocPt, NASSOCBINS, pT);

            if (iTriggBin >= 0) {
                histos->hPi0MassTrigg[iTriggBin]->Fill(mass);
            } else {
                for (int it = 0; it < NTRIGGBINS; it++) {
                    if (triggPt[it] < assocPt[iAssocBin+1]) continue;
                    if (binsWithTriggPeak[it] > 0 && iAssocBin >= 0) histos->hPi0MassAssocPeak[it][iAssocBin]->Fill(mass);
                    if (binsWithTriggSide[it] > 0 && iAssocBin >= 0) histos->hPi0MassAssocSide[it][iAssocBin]->Fill(mass);
                }
            }
        }
    }
}

void AliJHMRCorr::FillPionMassesLeading(TClonesArray *arrPhoton, TClonesArray *arrPi0, std::vector<int> listTrigg, detector idet, bool isPeakTriggLarger)
{
    if (listTrigg.size()<1) return;
    AliJBaseTrack *lvTrigg = (AliJBaseTrack*)arrPi0->At(listTrigg[0]);

    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)arrPhoton->At(i);
        if (!fIsFullSim && IsPhotonRemoved(lv1->E())) continue;
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)arrPhoton->At(j);
            if (GetAsymmetry(arrPhoton, lv1, lv2)>asymcut) continue;
            if (!fIsFullSim && IsPhotonRemoved(lv2->E())) continue;
            AliJBaseTrack lvSum = GetPhotonSumVector(arrPhoton, lv1, lv2);
            if (lvSum.Eta()<detEta[idet][0]+etacut || lvSum.Eta()>detEta[idet][1]-etacut) continue;
            double mass = 1000.*lvSum.M();
            double pT = lvSum.Pt();
            int iTriggBin = GetBin(triggPt, NTRIGGBINS, pT);
            int iAssocBin = GetBin(assocPt, NASSOCBINS, pT);
            if (iTriggBin >= 0 && pT==lvTrigg->Pt() && mass==1000.*lvTrigg->M()) {
                histos->hPi0MassTrigg[iTriggBin]->Fill(mass);
            } else {
                if (iAssocBin >= 0 && isPeakTriggLarger) histos->hPi0MassAssocPeak[iTriggBin][iAssocBin]->Fill(mass);
                if (iAssocBin >= 0 && !isPeakTriggLarger) histos->hPi0MassAssocSide[iTriggBin][iAssocBin]->Fill(mass);
            }
        }
    }
}

// Fill true pi0 mass for check
void AliJHMRCorr::FillPionMassesTrue(TClonesArray *arrPi0, int binsWithTriggReal[NTRIGGBINS], detector idet)
{
    int nPi0 = arrPi0->GetEntriesFast();
    for (int i = 1; i < nPi0; i++) {
        AliJBaseTrack *lvPi0 = (AliJBaseTrack*)arrPi0->At(i);
        if (lvPi0->Eta()<detEta[idet][0]+etacut || lvPi0->Eta()>detEta[idet][1]-etacut) continue;
        double mass = 1000.*lvPi0->M();
        double pT = lvPi0->Pt();
        int iTriggBin = GetBin(triggPt, NTRIGGBINS, pT);
        int iAssocBin = GetBin(assocPt, NASSOCBINS, pT);
        if (iTriggBin >= 0) histos->hPi0MassTrigg[iTriggBin]->Fill(mass);
        for (int it = 0; it < NTRIGGBINS; it++) {
            if (triggPt[it] < assocPt[iAssocBin+1]) continue;
            if (binsWithTriggReal[it] > 0 && iAssocBin >= 0) histos->hPi0MassAssocPeak[it][iAssocBin]->Fill(mass);
        }
    }
}

void AliJHMRCorr::FillAsymmetry(TClonesArray *arrPhoton, detector idet)
{
    int nPhoton = arrPhoton->GetEntriesFast();
    for (int i = 1; i < nPhoton; i++) {
        AliJBaseTrack *lv1 = (AliJBaseTrack*)arrPhoton->At(i);
        for (int j = 0; j < i; j++) {
            AliJBaseTrack *lv2 = (AliJBaseTrack*)arrPhoton->At(j);
            AliJBaseTrack lvSum = GetPhotonSumVector(arrPhoton, lv1, lv2);
            if (lvSum.Eta()<detEta[idet][0]+etacut || lvSum.Eta()>detEta[idet][1]-etacut) continue;
            double mass = 1000.*lvSum.M();
            if (!IsMassWindow(mass)) continue;
            double asym = GetAsymmetry(arrPhoton, lv1, lv2);
            histos->hEnergyAsymRec->Fill(asym);
            if (lv1->GetLabel()==1 && lv2->GetLabel()==1
                && lv1->GetMotherID()==lv2->GetMotherID())
                histos->hEnergyAsymTrue->Fill(asym);
        }
    }
}

// Return true if trigger and associated pi0s do not have common photons,
// false if either of the associated pi0 photons is used to reconstruct the
// trigger
bool AliJHMRCorr::CheckAssocPhotonPair(int iTrigg, int iAssoc, bool bMassWindow)
{
    std::vector<int> triggPairs;
    std::vector<int> assocPairs;
    if (bMassWindow) {
        if (photonId.size()<1) return true;
        triggPairs = photonId[iTrigg];
        assocPairs = photonId[iAssoc];
    } else {
        if (sidebandId.size()<1) return true;
        triggPairs = sidebandId[iTrigg];
        assocPairs = sidebandId[iAssoc];
    }
    for (int i = 0; i < 2; i++) {
        int triggId = triggPairs[i];
        for (int j = 0; j < 2; j++) {
            int assocId = assocPairs[j];
            if (triggId==assocId) return false;
        }
    }
    return true;
}
