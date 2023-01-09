#include "AnalyseCorrelations.h"

void AnalyseCorrelations()
{
    gStyle->SetOptStat(0);
    //processDataSTAR();
    processDataFoCal();
}

void processDataSTAR()
{
	TString fInName[ndata_star] = {
        "/home/heimarry/Simulations/focal/STAR-data/output_pAu.root",
        "/home/heimarry/Simulations/focal/STAR-data/output_pp.root"
        //"/home/heimarry/Simulations/focal-pythia-sim/star-pp_WRONG-TRUE-CORR-COMP.root"
	};

	TString fOutName[ndata_star] = {
        "analysis_STAR_pAu.root",
		"analysis_STAR_pp.root"
	};

    TString dataname[ndata_star] = {
        "STAR_pAu",
        "STAR_pp"
    };

	for (int idata = 0; idata < ndata_star; idata++) {
		fIn = TFile::Open(fInName[idata]);
		fOut = TFile::Open(fOutName[idata], "RECREATE");

		LoadInput();
        FitMassPeaks();
		GetScaleFactorsVersion1();
		DoAnalysis();
        DrawMassHistos(dataname[idata]);
		fOut->Write();

		fIn->Close();
		fOut->Close();
	}
}

void processDataFoCal()
{
	TString fInName[ndata_focal] = {
        //"/home/heimarry/Simulations/focal/analysis_output/2022-12-05_pp-focal_34.root"
        //"/home/heimarry/Simulations/focal/analysis_output/2022-12-14_pp-focal_no-pair-check.root"
        //"/home/heimarry/Simulations/focal/analysis_output/2022-12-16_pp-focal_pair-check_two-sidebands.root"
        //"/home/heimarry/Simulations/focal/analysis_output/2022-12-16_pp-focal_pair-check_two-sidebands.root"
        //"/home/heimarry/Simulations/focal/analysis_output/2022-12-19_pp-focal_no-pair-check.root"
        "/home/heimarry/Simulations/focal/analysis_output/2022-12-20_pp-focal_sideband-50-115-160-200.root"
	};

	TString fOutName[ndata_focal] = {
		"analysis_FoCal_pp.root"
	};

    TString dataname[ndata_focal] = {
        "FoCal_pp_fullsim"
    };

	for (int idata = 0; idata < ndata_focal; idata++) {
		fIn = TFile::Open(fInName[idata]);
		fOut = TFile::Open(fOutName[idata], "RECREATE");

		LoadInput();
        FitMassPeaks();
        //MixedEventCorrection();
		GetScaleFactorsVersion1();
		DoAnalysis();
        DrawMassHistos(dataname[idata]);
		fOut->Write();

		fIn->Close();
		fOut->Close();
	}
}

void LoadInput()
{

	hCounter = (TH1D*)fIn->Get("hCounter");
	nEvent = hCounter->GetBinContent(1);
	std::cout << "Input file contains " << nEvent << " events, proceed to analyse" << std::endl;

	hRealTriggCounter = (TH1D*)fIn->Get("hRealTriggCounter");
	for (int it = 0; it < nTriggBins; it++) nRealTrigg[it] = hRealTriggCounter->GetBinContent(it+1);

	for (int it = 0; it < nTriggBins; it++) {
	    double tlow = triggPt[it];
	    double tupp = triggPt[it+1];

	    hMassTrigg[it] = (TH1D*)fIn->Get(Form("Masses/hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp)); hMassTrigg[it]->Rebin(4);
        //hMassTrigg[it]->Scale(1., "width");

	    for (int ia = 0; ia < nAssocBins; ia++) {
	        double alow = assocPt[ia];
	        double aupp = assocPt[ia+1];

	        if (!useLeading && tlow < aupp) continue;

	        hMassAssocPeak[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hMassAssocPeak[it][ia]->Rebin(4);
            //hMassAssocPeak[it][ia]->Scale(1., "width");
	        hMassAssocSide[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hMassAssocSide[it][ia]->Rebin(4);
            //hMassAssocSide[it][ia]->Scale(1., "width");
            hMassAssocSum[it][ia] = (TH1D*)hMassAssocPeak[it][ia]->Clone(Form("hMassAssocSum[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hMassAssocSum[it][ia]->Rebin(4);
            hMassAssocSum[it][ia]->Add(hMassAssocSide[it][ia]);

	        hCorrReal[it][ia]     = (TH2D*)fIn->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));           //hCorrReal[it][ia]->Rebin2D(12);
	        hCorrMassMass[it][ia] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrMassMass[it][ia]->Rebin2D(12);
            hCorrMassSide[it][ia] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrMassSide[it][ia]->Rebin2D(12);
	        hCorrSideMass[it][ia] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrSideMass[it][ia]->Rebin2D(12);
	        hCorrSideSide[it][ia] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrSideSide[it][ia]->Rebin2D(12);
            hCorrMeas[it][ia]     = (TH2D*)hCorrMassMass[it][ia]->Clone(Form("hCorrMeas[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

            hCorrMassMassMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrMassMassMixed[it][ia]->Rebin2D(12);
            hCorrMassSideMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrMassSideMixed[it][ia]->Rebin2D(12);
            hCorrSideMassMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrSideMassMixed[it][ia]->Rebin2D(12);
            hCorrSideSideMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hCorrSideSideMixed[it][ia]->Rebin2D(12);
            hCorrMeasMixed[it][ia]     = (TH2D*)hCorrMassMassMixed[it][ia]->Clone(Form("hCorrMeasMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

	    }
	}
}

void DoAnalysis()
{
    for (int it = 0; it < nTriggBins; it++) {
        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (!useLeading && tlow < aupp) continue;

            hCorrRealProj[it][ia] = hCorrReal[it][ia]->ProjectionX();
            //hCorrRealProj[it][ia]->Scale(1.0/nRealTrigg[it], "width");
            //hCorrRealProj[it][ia]->Scale(1.0/hCorrRealProj[it][ia]->GetEntries(), "width");

            hCorrMassMassProj[it][ia] = hCorrMassMass[it][ia]->ProjectionX();

            hCorrMassSideProj[it][ia] = hCorrMassSide[it][ia]->ProjectionX();
            hCorrMassSideProj[it][ia]->Scale(alpha[it][ia]);

            hCorrSideMassProj[it][ia] = hCorrSideMass[it][ia]->ProjectionX();
            hCorrSideMassProj[it][ia]->Scale(beeta[it][ia]);

            hCorrSideSideProj[it][ia] = hCorrSideSide[it][ia]->ProjectionX();
            hCorrSideSideProj[it][ia]->Scale(yamma[it][ia]);

            hCorr[it][ia] = (TH1D*)hCorrMassMassProj[it][ia]->Clone(Form("hCorrFinal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorr[it][ia]->Add(hCorrMassSideProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideMassProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideSideProj[it][ia]);

            // Take efficiencies into account
            //hCorr[it][ia]->Scale(1.0/st[it], "width");
            //hCorr[it][ia]->Scale(1.0/hCorr[it][ia]->GetEntries(), "width");

            double massWindowMin, massWindowMax;
            if (bUseConstMassWindow) {
                massWindowMin = massMin;
                massWindowMax = massMax;
            } else {
                massWindowMin = massPeakPosTrigg[it]-3.*massSigmaTrigg[it];
                massWindowMax = massPeakPosTrigg[it]+3.*massSigmaTrigg[it];
            }
            hCorrMeasProj[it][ia] = hCorrMeas[it][ia]->ProjectionX();
            //hCorrMeasProj[it][ia]->Scale(1.0/fFitTrigg[it]->Integral(massWindowMin, massWindowMax), "width");
            //hCorrMeasProj[it][ia]->Scale(1.0/hMassTrigg[it]->Integral(hMassTrigg[it]->GetXaxis()->FindBin(massWindowMin), hMassTrigg[it]->GetXaxis()->FindBin(massWindowMax)), "width");
            //hCorrMeasProj[it][ia]->Scale(1.0/hCorrMeasProj[it][ia]->GetEntries(), "width");
        }
    }

}

void MixedEventCorrection()
{
    double alpha = 0.;
    for (int it = 0; it < nTriggBins; it++) {
        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (!useLeading && tlow < aupp) continue;

            alpha = GetMixedEventNormalization(hCorrMeasMixed[it][ia]);
            hCorrMeasMixed[it][ia]->Scale(alpha);
            hCorrMeas[it][ia]->Divide(hCorrMeasMixed[it][ia]);

            alpha = GetMixedEventNormalization(hCorrMassMassMixed[it][ia]);
            hCorrMassMassMixed[it][ia]->Scale(alpha);
            hCorrMassMass[it][ia]->Divide(hCorrMassMassMixed[it][ia]);

            alpha = GetMixedEventNormalization(hCorrMassSideMixed[it][ia]);
            hCorrMassSideMixed[it][ia]->Scale(alpha);
            hCorrMassSide[it][ia]->Divide(hCorrMassSideMixed[it][ia]);

            alpha = GetMixedEventNormalization(hCorrSideMassMixed[it][ia]);
            hCorrSideMassMixed[it][ia]->Scale(alpha);
            hCorrSideMass[it][ia]->Divide(hCorrSideMassMixed[it][ia]);

            alpha = GetMixedEventNormalization(hCorrSideSideMixed[it][ia]);
            hCorrSideSideMixed[it][ia]->Scale(alpha);
            hCorrSideSide[it][ia]->Divide(hCorrSideSideMixed[it][ia]);
        }
    }
}

double GetMixedEventNormalization(TH2D* h)
{
    double norm = 0;
    int binx = h->GetXaxis()->FindBin(0.);
    int biny = h->GetYaxis()->FindBin(0.);
    norm += h->GetBinContent(binx, biny);
    norm += h->GetBinContent(binx+1, biny);
    norm += h->GetBinContent(binx-1, biny);
    norm += h->GetBinContent(binx, biny+1);
    norm += h->GetBinContent(binx, biny-1);
    norm /= 5.;
    return 1./norm;
}

void FitMassPeaks()
{
    for (int it = 0; it < nTriggBins; it++) {
        fFitTrigg[it] = new TF1(Form("fFitTrigg%d", it), FitFunction, 10, 500, 10);
        fFitTrigg[it]->SetParameters(0., 0., 0., 0., 0., 1., 1., 135., 5., 10.);
        //fFitTrigg[it]->SetParLimits(4, 132., 137.);
        fFitTrigg[it]->SetParLimits(7, 130., 160.);
        fFitTrigg[it]->SetParLimits(8, 0.1, 50.);
        fFitTrigg[it]->SetParLimits(9, 5., 50.);
        fFitTrigg[it]->FixParameter(0, 0.);
        fFitTrigg[it]->FixParameter(1, 0.);
        fFitTrigg[it]->FixParameter(6, 0.);
        fFitTrigg[it]->FixParameter(9, 0.);
        fFitTrigg[it]->SetNpx(1000);

        fPeakTrigg[it] = new TF1(Form("fPeakTrigg%d", it), FitPeak, 10, 500, 5);
        fPeakTrigg[it]->SetLineColor(kBlue);
        fPeakTrigg[it]->SetNpx(1000);

        fBgTrigg[it] = new TF1(Form("fBgTrigg%d", it), FitBackground, 10, 500, 5);
        fBgTrigg[it]->SetLineColor(kBlack);
        fBgTrigg[it]->SetLineStyle(kDashed);
        fBgTrigg[it]->SetNpx(1000);

        //hMassTrigg[it]->Fit(Form("fFitTrigg%d", it), "QNRW+");
        hMassTrigg[it]->Fit(Form("fFitTrigg%d", it), "SQNR+");
        fitConstantVal[it] = fFitTrigg[it]->GetParameter(3);
        fFitTrigg[it]->GetParameters(parTrigg);
        fBgTrigg[it]->SetParameters(parTrigg);
        fPeakTrigg[it]->SetParameters(&parTrigg[5]);

        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (!useLeading && tlow < aupp) continue;

            // Assoc mass histrograms where trigger is from peak region
            fFitAssocPeak[it][ia] = new TF1(Form("fFitAssocPeak%d-%d", it, ia), FitFunction, 50, 500, 10);
            fFitAssocPeak[it][ia]->SetParameters(0., 0., 0., 0., 0., 1., 1., 135., 5., 10.);
            //fFitAssocPeak[it]->SetParLimits(4, 132., 137.);
            fFitAssocPeak[it][ia]->SetParLimits(7, 130., 160.);
            fFitAssocPeak[it][ia]->SetParLimits(8, 0.1, 50.);
            fFitAssocPeak[it][ia]->SetParLimits(9, 5., 50.);
            fFitAssocPeak[it][ia]->FixParameter(0, 0.);
            fFitAssocPeak[it][ia]->FixParameter(1, 0.);
            fFitAssocPeak[it][ia]->FixParameter(6, 0.);
            fFitAssocPeak[it][ia]->FixParameter(9, 0.);
            fFitAssocPeak[it][ia]->SetNpx(1000);

            fPeakAssocPeak[it][ia] = new TF1(Form("fPeakAssocPeak%d-%d", it, ia), FitPeak, 20, 500, 5);
            fPeakAssocPeak[it][ia]->SetLineColor(kBlue);
            fPeakAssocPeak[it][ia]->SetNpx(1000);

            fBgAssocPeak[it][ia]= new TF1(Form("fBgAssocPeak%d-%d", it, ia), FitBackground, 20, 500, 5);
            fBgAssocPeak[it][ia]->SetLineColor(kBlack);
            fBgAssocPeak[it][ia]->SetLineStyle(kDashed);
            fBgAssocPeak[it][ia]->SetNpx(1000);

            hMassAssocPeak[it][ia]->SetOption("p");
            hMassAssocPeak[it][ia]->Fit(Form("fFitAssocPeak%d-%d", it, ia), "SQNR+");
            fFitAssocPeak[it][ia]->GetParameters(parAssocPeak);
            fBgAssocPeak[it][ia]->SetParameters(parAssocPeak);
            fPeakAssocPeak[it][ia]->SetParameters(&parAssocPeak[5]);

            // Assoc mass histrograms where trigger is from sideband region
            fFitAssocSide[it][ia] = new TF1(Form("fFitAssocSide%d-%d", it, ia), FitFunction, 50, 500, 10);
            fFitAssocSide[it][ia]->SetParameters(0., 0., 0., 0., 0., 1., 1., 135., 5., 10.);
            //fFitAssocSide[it]->SetParLimits(4, 132., 137.);
            fFitAssocSide[it][ia]->SetParLimits(7, 130., 160.);
            fFitAssocSide[it][ia]->SetParLimits(8, 0.1, 50.);
            fFitAssocSide[it][ia]->SetParLimits(9, 5., 50.);
            fFitAssocSide[it][ia]->FixParameter(0, 0.);
            fFitAssocSide[it][ia]->FixParameter(1, 0.);
            fFitAssocSide[it][ia]->FixParameter(6, 0.);
            fFitAssocSide[it][ia]->FixParameter(9, 0.);
            fFitAssocSide[it][ia]->SetNpx(1000);

            fPeakAssocSide[it][ia] = new TF1(Form("fPeakAssocSide%d-%d", it, ia), FitPeak, 50, 500, 5);
            fPeakAssocSide[it][ia]->SetLineColor(kBlue);
            fPeakAssocSide[it][ia]->SetNpx(1000);

            fBgAssocSide[it][ia]= new TF1(Form("fBgAssocSide%d-%d", it, ia), FitBackground, 50, 500, 5);
            fBgAssocSide[it][ia]->SetLineColor(kBlack);
            fBgAssocSide[it][ia]->SetLineStyle(kDashed);
            fBgAssocSide[it][ia]->SetNpx(1000);

            hMassAssocSide[it][ia]->Fit(Form("fFitAssocSide%d-%d", it, ia), "SQNR+");
            fFitAssocSide[it][ia]->GetParameters(parAssocSide);
            fBgAssocSide[it][ia]->SetParameters(parAssocSide);
            fPeakAssocSide[it][ia]->SetParameters(&parAssocSide[5]);

            // Assoc mass histrograms where trigger is from their sum
            fFitAssocSum[it][ia] = new TF1(Form("fFitAssocSum%d-%d", it, ia), FitFunction, 50, 500, 10);
            fFitAssocSum[it][ia]->SetParameters(0., 0., 0., 0., 0., 1., 1., 135., 5., 10.);
            fFitAssocSum[it][ia]->SetParLimits(7, 130., 160.);
            fFitAssocSum[it][ia]->SetParLimits(8, 0.1, 50.);
            fFitAssocSum[it][ia]->SetParLimits(9, 5., 50.);
            fFitAssocSum[it][ia]->FixParameter(0, 0.);
            fFitAssocSum[it][ia]->FixParameter(1, 0.);
            fFitAssocSum[it][ia]->FixParameter(6, 0.);
            fFitAssocSum[it][ia]->FixParameter(9, 0.);
            fFitAssocSum[it][ia]->SetNpx(1000);

            fPeakAssocSum[it][ia] = new TF1(Form("fPeakAssocSum%d-%d", it, ia), FitPeak, 50, 500, 5);
            fPeakAssocSum[it][ia]->SetLineColor(kBlue);
            fPeakAssocSum[it][ia]->SetNpx(1000);

            fBgAssocSum[it][ia]= new TF1(Form("fBgAssocSum%d-%d", it, ia), FitBackground, 50, 500, 5);
            fBgAssocSum[it][ia]->SetLineColor(kBlack);
            fBgAssocSum[it][ia]->SetLineStyle(kDashed);
            fBgAssocSum[it][ia]->SetNpx(1000);

            hMassAssocSum[it][ia]->Fit(Form("fFitAssocSum%d-%d", it, ia), "SQNR+");
            fFitAssocSum[it][ia]->GetParameters(parAssocSum);
            fBgAssocSum[it][ia]->SetParameters(parAssocSum);
            fPeakAssocSum[it][ia]->SetParameters(&parAssocSum[5]);
        }
    }
}

void GetScaleFactorsVersion1()
{
    double massWindowMin, massWindowMax;
    for (int it = 0; it < nTriggBins; it++) {
        if (bUseConstMassWindow) {
            massWindowMin = massMin;
            massWindowMax = massMax;
        } else {
            massWindowMin = massPeakPosTrigg[it]-3.*massSigmaTrigg[it];
            massWindowMax = massPeakPosTrigg[it]+3.*massSigmaTrigg[it];
        }
        //st[it] = fPeakTrigg[it]->Integral(massWindowMin, massWindowMax);
        st[it]  = hMassTrigg[it]->Integral(hMassTrigg[it]->GetXaxis()->FindBin(massWindowMin), hMassTrigg[it]->GetXaxis()->FindBin(massWindowMax)) - fBgTrigg[it]->Integral(massWindowMin, massWindowMax)/hMassTrigg[it]->GetBinWidth(0);
        //st[it] *= effCorrTrigg[it]*pi0br;
        stobtrigg[it] = fPeakTrigg[it]->Integral(massWindowMin, massWindowMax)/fFitTrigg[it]->Integral(massWindowMin, massWindowMax);
        std::cout << "\n\tbin [ " << triggPt[it] << " " << triggPt[it+1] << " ] : \treal=" << nRealTrigg[it] << "\treconst=" << st[it] << "\trec=" << stfit[it] << "\trec/real=" << st[it]/nRealTrigg[it] << std::endl;
        std::cout << "\tS/(S+B) : " << stobtrigg[it] << std::endl;
        //std::cout << "bg : " << fBgTrigg[it]->Integral(massWindowMin, massWindowMax) << std::endl;
        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];
            if (bUseConstMassWindow) {
                massWindowMin = massMin;
                massWindowMax = massMax;
            } else {
                massWindowMin = massPeakPosAssoc[it]-3.*massSigmaAssoc[it];
                massWindowMax = massPeakPosAssoc[it]+3.*massSigmaAssoc[it];
            }

            if (!useLeading && tlow < aupp) continue;
            //alpha[it][ia]  = fBgAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax)/(fBgAssocPeak[it][ia]->Integral(40, 80) + fBgAssocPeak[it][ia]->Integral(210, 280));
            //beeta[it][ia]  = fBgTrigg[it]->Integral(massWindowMin, massWindowMax)/(fBgTrigg[it]->Integral(40, 80) + fBgTrigg[it]->Integral(210, 280));

            alpha[it][ia]  = fBgAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax)/(fBgAssocPeak[it][ia]->Integral(50, 115) + fBgAssocPeak[it][ia]->Integral(160, 200));
            beeta[it][ia]  = fBgTrigg[it]->Integral(massWindowMin, massWindowMax)/(fBgTrigg[it]->Integral(50, 115) + fBgTrigg[it]->Integral(160, 200));

            //alpha[it][ia] = fBgAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax)/fBgAssocPeak[it][ia]->Integral(300, 450);
            //beeta[it][ia] = fBgTrigg[it]->Integral(massWindowMin, massWindowMax)/fBgTrigg[it]->Integral(300, 450);
            //cout << fBgAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax) << "/" << fBgAssocPeak[it][ia]->Integral(300, 450) << endl;
            yamma[it][ia] = alpha[it][ia]*beeta[it][ia];
            std::cout << "\t\tbin [ "  << assocPt[ia] << " " << assocPt[ia+1] << " ] : "
                      << "\talpha=" << alpha[it][ia]
                      << "\tbeeta="  << beeta[it][ia]
                      << "\tgamma=" << yamma[it][ia] << std::endl;
            stobassoc[it][ia] = fPeakAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax)/fFitAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax);
            std::cout << "\t\tS/(S+B) : " << stobassoc[it][ia] << std::endl;
        }
    }
}

void GetScaleFactorsVersion2()
{
    double massWindowMin, massWindowMax;
    for (int it = 0; it < nTriggBins; it++) {
        //st[it] = fPeakTrigg[it]->Integral(massWindowMin, massWindowMax);
        if (bUseConstMassWindow) {
            massWindowMin = massMin;
            massWindowMax = massMax;
        } else {
            massWindowMin = massPeakPosTrigg[it]-3.*massSigmaTrigg[it];
            massWindowMax = massPeakPosTrigg[it]+3.*massSigmaTrigg[it];
        }
        //st[it] = fPeakTrigg[it]->Integral(massWindowMin, massWindowMax);
        st[it]  = hMassTrigg[it]->Integral(hMassTrigg[it]->GetXaxis()->FindBin(massWindowMin), hMassTrigg[it]->GetXaxis()->FindBin(massWindowMax)) - fBgTrigg[it]->Integral(massWindowMin, massWindowMax)/hMassTrigg[it]->GetBinWidth(0);
        std::cout << "\n\tbin [ " << triggPt[it] << " " << triggPt[it+1] << " ] : \treal=" << nRealTrigg[it] << "\treconst=" << st[it] << "\trec/real=" << st[it]/nRealTrigg[it] << std::endl;
        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];
            if (bUseConstMassWindow) {
                massWindowMin = massMin;
                massWindowMax = massMax;
            } else {
                massWindowMin = massPeakPosAssoc[it]-3.*massSigmaAssoc[it];
                massWindowMax = massPeakPosAssoc[it]+3.*massSigmaAssoc[it];
            }

            if (!useLeading && tlow < aupp) continue;
            double alpha1 = fBgAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax)/fBgAssocPeak[it][ia]->Integral(300, 350);
            double alpha2 = fBgAssocSide[it][ia]->Integral(massWindowMin, massWindowMax)/fBgAssocSide[it][ia]->Integral(300, 350);
            double alpha3 = fBgAssocPeak[it][ia]->Integral(300, 350)/fBgAssocSide[it][ia]->Integral(300, 350);
            //double alpha1 = fBgAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax)/(fBgAssocPeak[it][ia]->Integral(40, 80) + fBgAssocPeak[it][ia]->Integral(210, 280));
            //double alpha2 = fBgAssocSide[it][ia]->Integral(massWindowMin, massWindowMax)/(fBgAssocSide[it][ia]->Integral(40, 80) + fBgAssocSide[it][ia]->Integral(210, 280));
            //double alpha3 = (fBgAssocPeak[it][ia]->Integral(40, 80) + fBgAssocPeak[it][ia]->Integral(210, 280))/(fBgAssocSide[it][ia]->Integral(40, 80) + fBgAssocSide[it][ia]->Integral(210, 280));
            double alpha4 = hMassAssocPeak[it][ia]->Integral(hMassAssocPeak[it][ia]->GetXaxis()->FindBin(massWindowMin), hMassAssocPeak[it][ia]->GetXaxis()->FindBin(massWindowMax)) - fBgAssocPeak[it][ia]->Integral(massWindowMin, massWindowMax)/hMassAssocPeak[it][ia]->GetBinWidth(0);
            alpha4 /= hMassAssocSide[it][ia]->Integral(hMassAssocSide[it][ia]->GetXaxis()->FindBin(massWindowMin), hMassAssocSide[it][ia]->GetXaxis()->FindBin(massWindowMax)) - fBgAssocSide[it][ia]->Integral(massWindowMin, massWindowMax)/hMassAssocSide[it][ia]->GetBinWidth(0);
            double beeta1 = fBgTrigg[it]->Integral(massWindowMin, massWindowMax)/fBgTrigg[it]->Integral(300, 450);
            //double beeta1 = fBgTrigg[it]->Integral(massWindowMin, massWindowMax)/(fBgTrigg[it]->Integral(40, 80) + fBgTrigg[it]->Integral(210, 280));

            alpha[it][ia] = alpha1;
            beeta[it][ia] = alpha4*beeta1;
            yamma[it][ia] = (alpha3 - alpha1*alpha3 - alpha2*alpha4)*beeta1;
            std::cout << "\t\tbin [ "  << assocPt[ia] << " " << assocPt[ia+1] << " ] : "
                      << "\talpha1=" << alpha1
                      << "\talpha2=" << alpha2
                      << "\talpha3=" << alpha3
                      << "\talpha4=" << alpha4
                      << "\tbeta1=" << beeta1
                      << "\talpha=" << alpha[it][ia]
                      << "\tbeeta="  << beeta[it][ia]
                      << "\tgamma=" << yamma[it][ia] << std::endl;
        }
    }
}

void DrawMassHistos(TString dataname)
{
    TLatex tl;
    tl.SetTextSize(0.035);
    tl.SetTextFont(42);
    double massWindowMin, massWindowMax;
    for (int it = 0; it < nTriggBins; it++) {
        if (bUseConstMassWindow) {
            massWindowMin = massMin;
            massWindowMax = massMax;
        } else {
            massWindowMin = massPeakPosTrigg[it]-3.*massSigmaTrigg[it];
            massWindowMax = massPeakPosTrigg[it]+3.*massSigmaTrigg[it];
        }
        cMassTrigg[it] = new TCanvas(Form("cMassTrigg%d", it), Form("cMassTrigg%d", it), 600, 600);
        hMassTrigg[it]->SetTitle(Form("Trigger invariant mass, %.1f < p_{T,t} < %.1f GeV/c; M_{#gamma#gamma} (MeV/c^{2}); counts", triggPt[it], triggPt[it+1]));
        hMassTrigg[it]->GetXaxis()->SetRangeUser(0., 500.);
        hMassTrigg[it]->GetYaxis()->SetMaxDigits(3);

        fPeakColored[it] = (TF1*)fFitTrigg[it]->Clone("fPeakColored");
        fPeakColored[it]->SetFillColor(kRed-10);
        fPeakColored[it]->GetXaxis()->SetRangeUser(massWindowMin, massWindowMax);
        fPeakColored[it]->SetFillStyle(1001);
        fLeftSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fLeftSidebandColored[it]->SetFillColor(kGreen-10);
        fLeftSidebandColored[it]->GetXaxis()->SetRangeUser(40., 80.);
        fLeftSidebandColored[it]->SetFillStyle(1001);
        fRightSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fRightSidebandColored[it]->SetFillColor(kGreen-10);
        fRightSidebandColored[it]->GetXaxis()->SetRangeUser(300., 450.);
        fRightSidebandColored[it]->SetFillStyle(1001);
        fBgColored[it] = (TF1*)fBgTrigg[it]->Clone("fBgColored");
        fBgColored[it]->SetFillColor(kBlue-10);
        fBgColored[it]->GetXaxis()->SetRangeUser(massWindowMin, massWindowMax);
        fBgColored[it]->SetFillStyle(1001);

        cMassTrigg[it]->cd();
        hMassTrigg[it]->Draw("HIST");
        fPeakColored[it]->Draw("SAME");
        //fLeftSidebandColored[it]->Draw("SAME");
        fRightSidebandColored[it]->Draw("SAME");
        fBgColored[it]->Draw("SAME");
        fFitTrigg[it]->Draw("SAME");
        fBgTrigg[it]->Draw("SAME");
        fPeakTrigg[it]->Draw("SAME");

        tl.DrawLatexNDC(.52, .7, "#pi^{0} mass peak");
        tl.DrawLatexNDC(.55, .65, Form("#bf{#bf{%.1f < p_{T,t} < %.1f GeV/c}}", triggPt[it], triggPt[it+1]));
        gPad->RedrawAxis();
        cMassTrigg[it]->SaveAs(Form("MassTrigg%d_%s.pdf", it, dataname.Data()));

        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];
            if (bUseConstMassWindow) {
                massWindowMin = massMin;
                massWindowMax = massMax;
            } else {
                massWindowMin = massPeakPosAssoc[it]-3.*massSigmaAssoc[it];
                massWindowMax = massPeakPosAssoc[it]+3.*massSigmaAssoc[it];
            }

            if (!useLeading && tlow < aupp) continue;

            // Trigger from mass region
            cMassAssocPeak[it][ia] = new TCanvas(Form("cMassAssocPeak%d:%d",it,ia), Form("cMassAssocPeak%d:%d",it,ia), 600, 600);
            hMassAssocPeak[it][ia]->SetTitle("; M_{#gamma#gamma} (MeV/c^{2}); counts");
            //hMassAssocPeak[it][ia]->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; M_{#gamma#gamma}; counts", tlow, tupp, alow, aupp));
            hMassAssocPeak[it][ia]->GetXaxis()->SetRangeUser(0., 500.);
            hMassAssocPeak[it][ia]->GetYaxis()->SetMaxDigits(3);

            fPeakColoredAssocPeak[it][ia] = (TF1*)fFitAssocPeak[it][ia]->Clone("fPeakColoredAssocPeak");
            fPeakColoredAssocPeak[it][ia]->SetFillColor(kRed-10);
            fPeakColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(massWindowMin, massWindowMax);
            fPeakColoredAssocPeak[it][ia]->SetFillStyle(1001);
            fLeftSidebandColoredAssocPeak[it][ia] = (TF1*)fFitAssocPeak[it][ia]->Clone("fLeftSidebandColoredAssocPeak");
            fLeftSidebandColoredAssocPeak[it][ia]->SetFillColor(kGreen-10);
            fLeftSidebandColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(40., 80.);
            fLeftSidebandColoredAssocPeak[it][ia]->SetFillStyle(1001);
            fRightSidebandColoredAssocPeak[it][ia] = (TF1*)fFitAssocPeak[it][ia]->Clone("fLeftSidebandColoredAssocPeak");
            fRightSidebandColoredAssocPeak[it][ia]->SetFillColor(kGreen-10);
            fRightSidebandColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(300., 450.);
            fRightSidebandColoredAssocPeak[it][ia]->SetFillStyle(1001);
            fBgColoredAssocPeak[it][ia] = (TF1*)fBgAssocPeak[it][ia]->Clone("fBgColoredAssocPeak");
            fBgColoredAssocPeak[it][ia]->SetFillColor(kBlue-10);
            fBgColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(massWindowMin, massWindowMax);
            fBgColoredAssocPeak[it][ia]->SetFillStyle(1001);

            cMassAssocPeak[it][ia]->cd();
            hMassAssocPeak[it][ia]->Draw("HIST");
            fPeakColoredAssocPeak[it][ia]->Draw("SAME");
            //fLeftSidebandColoredAssocPeak[it][ia]->Draw("SAME");
            fRightSidebandColoredAssocPeak[it][ia]->Draw("SAME");
            fBgColoredAssocPeak[it][ia]->Draw("SAME");
            hMassAssocPeak[it][ia]->Draw("HIST SAME");
            fFitAssocPeak[it][ia]->Draw("SAME");
            fBgAssocPeak[it][ia]->Draw("SAME");
            fPeakAssocPeak[it][ia]->Draw("SAME");

            tl.DrawLatexNDC(.52, .7, "#pi^{0} mass peak");
            tl.DrawLatexNDC(.55, .65, Form("#bf{#bf{%.1f < p_{T,t} < %.1f GeV/c}}", tlow, tupp));
            tl.DrawLatexNDC(.55, .6, Form("%.1f < p_{T,a} < %.1f GeV/c", alow, aupp));
            gPad->RedrawAxis();
            cMassAssocPeak[it][ia]->SaveAs(Form("MassAssocPeak%d-%d_%s.pdf",it,ia,dataname.Data()));

            // Trigger from sideband region
            cMassAssocSide[it][ia] = new TCanvas(Form("cMassAssocSide%d:%d",it,ia), Form("cMassAssocSide%d:%d",it,ia), 600, 600);
            //hMassAssocSide[it][ia]->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; M_{#gamma#gamma}; counts", tlow, tupp, alow, aupp));
            hMassAssocSide[it][ia]->SetTitle("; M_{#gamma#gamma}; counts");
            hMassAssocSide[it][ia]->GetXaxis()->SetRangeUser(0., 500.);
            hMassAssocSide[it][ia]->GetYaxis()->SetMaxDigits(6);

            fPeakColoredAssocSide[it][ia] = (TF1*)fFitAssocSide[it][ia]->Clone("fPeakColoredAssocSide");
            fPeakColoredAssocSide[it][ia]->SetFillColor(kRed-10);
            fPeakColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(massWindowMin, massWindowMax);
            fPeakColoredAssocSide[it][ia]->SetFillStyle(1001);
            fLeftSidebandColoredAssocSide[it][ia] = (TF1*)fFitAssocSide[it][ia]->Clone("fLeftSidebandColoredAssocSide");
            fLeftSidebandColoredAssocSide[it][ia]->SetFillColor(kGreen-10);
            fLeftSidebandColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(40., 80.);
            fLeftSidebandColoredAssocSide[it][ia]->SetFillStyle(1001);
            fRightSidebandColoredAssocSide[it][ia] = (TF1*)fFitAssocSide[it][ia]->Clone("fLeftSidebandColoredAssocSide");
            fRightSidebandColoredAssocSide[it][ia]->SetFillColor(kGreen-10);
            fRightSidebandColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(300., 450.);
            fRightSidebandColoredAssocSide[it][ia]->SetFillStyle(1001);
            fBgColoredAssocSide[it][ia] = (TF1*)fBgAssocSide[it][ia]->Clone("fBgColoredAssocSide");
            fBgColoredAssocSide[it][ia]->SetFillColor(kBlue-10);
            fBgColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(massWindowMin, massWindowMax);
            fBgColoredAssocSide[it][ia]->SetFillStyle(1001);

            cMassAssocSide[it][ia]->cd();
            hMassAssocSide[it][ia]->Draw("HIST");
            fPeakColoredAssocSide[it][ia]->Draw("SAME");
            //fLeftSidebandColoredAssocSide[it][ia]->Draw("SAME");
            fRightSidebandColoredAssocSide[it][ia]->Draw("SAME");
            fBgColoredAssocSide[it][ia]->Draw("SAME");
            fFitAssocSide[it][ia]->Draw("SAME");
            fBgAssocSide[it][ia]->Draw("SAME");
            fPeakAssocSide[it][ia]->Draw("SAME");

            tl.DrawLatexNDC(.52, .7, "#pi^{0} mass peak");
            tl.DrawLatexNDC(.55, .65, Form("#bf{#bf{%.1f < p_{T,t} < %.1f GeV/c}}", tlow, tupp));
            tl.DrawLatexNDC(.55, .6, Form("%.1f < p_{T,a} < %.1f GeV/c", alow, aupp));
            gPad->RedrawAxis();
            //cMassAssocSide[it][ia]->SaveAs(Form("MassAssocSide%d-%d_%s.pdf",it,ia,dataname.Data()));
        }
    }
}

double FitPeak(double *x, double *p)
{
    double c1 = p[0];
    double c2 = p[1];
    double mu = p[2];
    double sigma1 = p[3];
    double sigma2 = p[4];

    return c1*TMath::Gaus(x[0], mu, sigma1) + c2*TMath::Gaus(x[0], mu, sigma2);

    //return c1*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma1*sigma1)) + c2*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma2*sigma2));
}

double FitBackground(double *x, double *p)
{
    double a = p[0];
    double b = p[1];
    double c = p[2];
    double d = p[3];
    double e = p[4];

    return a*x[0]*x[0]*x[0]*x[0] + b*x[0]*x[0]*x[0] + c*x[0]*x[0] + d*x[0] + e;
}

double FitFunction(double *x, double *p)
{
    return FitBackground(x, p) + FitPeak(x, &p[5]);
}
