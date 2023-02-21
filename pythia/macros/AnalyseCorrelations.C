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
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-02-09_pp-focal_smaller-mw_no-asym.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-02-09_pp-focal_smaller-mw_asym-08.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-02-09_pp-focal_smaller-mw_asym-05.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-02-09_pp-focal_smaller-mw_asym-08_w.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-02-10_pp-focal_smaller-mw_asym-05_w.root"
        //************************************************************************************************
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-01-26_pp-focal_no-photon-eff.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-01-27_pp-focal_asym-08.root"
        "/home/heimarry/Simulations/focal/analysis_output_2/2023-02-14_gAnalysis_away-side_2.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-01-30_gAnalysis_asym-08_etacut-02_weighted.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-01-27_pp-focal_asym-08_weight-078.root"
        //"/home/heimarry/Simulations/focal/analysis_output_2/2023-02-01_pp-focal_asym-05_weight-05_2.root"
        //"/home/heimarry/Simulations/focal-pythia-sim/focal-pp_no-weight.root"
	};

	TString fOutName[ndata_focal] = {
		"analysis_FoCal_pp.root"
		//"analysis_FoCal_pp_full-sim.root"
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
        //DoSimplifiedAnalysis();
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

	    hMassTrigg[it] = (TH1D*)fIn->Get(Form("Masses/hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp));
        //hMassTrigg[it]->Rebin(5);

	    for (int ia = 0; ia < nAssocBins; ia++) {
	        double alow = assocPt[ia];
	        double aupp = assocPt[ia+1];

	        if (tlow < aupp) continue;

	        hMassAssocPeak[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            //hMassAssocPeak[it][ia]->Rebin(5);
	        hMassAssocSide[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            //hMassAssocSide[it][ia]->Rebin(5);

            // 2D same event correlations
	        hCorrReal[it][ia]     = (TH2D*)fIn->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));           
	        hCorrMassMass[it][ia] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassSide[it][ia] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
	        hCorrSideMass[it][ia] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
	        hCorrSideSide[it][ia] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

            // Prjections of correlation functions
            hCorrRealProj[it][ia] = hCorrReal[it][ia]->ProjectionX();
            hCorrMassMassProj[it][ia] = hCorrMassMass[it][ia]->ProjectionX();
            hCorrMassSideProj[it][ia] = hCorrMassSide[it][ia]->ProjectionX();
            hCorrSideMassProj[it][ia] = hCorrSideMass[it][ia]->ProjectionX();
            hCorrSideSideProj[it][ia] = hCorrSideSide[it][ia]->ProjectionX();
            hCorrNonCorrected[it][ia] = (TH1D*)hCorrMassMassProj[it][ia]->Clone(Form("hCorrNonCorrected[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrCorrected[it][ia] = (TH1D*)hCorrMassMassProj[it][ia]->Clone(Form("hCorrCorrected[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

            // Rebin if needed
            hCorrRealProj[it][ia]->Rebin(8);
            hCorrMassMassProj[it][ia]->Rebin(8);
            hCorrMassSideProj[it][ia]->Rebin(8);
            hCorrSideMassProj[it][ia]->Rebin(8);
            hCorrSideSideProj[it][ia]->Rebin(8);
            hCorrNonCorrected[it][ia]->Rebin(8);
            hCorrCorrected[it][ia]->Rebin(8);

            // 2D mixed event correlations
            hCorrMassMassMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassSideMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideMassMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideSideMixed[it][ia] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
	    }
	}
}

void DoAnalysis()
{
    for (int it = 0; it < nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];
        int nreal = nRealTrigg[it];
        int ncandidates = hMassTrigg[it]->Integral(hMassTrigg[it]->GetXaxis()->FindBin(massMin), hMassTrigg[it]->GetXaxis()->FindBin(massMax));
        int ncandidateswobg = st[it];
        cout << "Trigg bin [" << tlow << " " << tupp << "] : \tN REAL TRIGGERS=" << nreal << "\tN CANDIDATES="  << ncandidates << "\tN CANDITATES (NO BG)=" << ncandidateswobg << endl;
        for (int ia = 0; ia < nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            // True correlation
            hCorrRealProj[it][ia]->Scale(1.0/nreal, "width");

            // Non sideband corrected correlation
            hCorrNonCorrected[it][ia]->Scale(1.0/ncandidates, "width");

            // Sideband corrected correlation
            hCorrMassSideProj[it][ia]->Scale(alpha[it][ia]);
            hCorrSideMassProj[it][ia]->Scale(beeta[it][ia]);
            hCorrSideSideProj[it][ia]->Scale(yamma[it][ia]);

            hCorrCorrected[it][ia]->Add(hCorrMassSideProj[it][ia], -1);
            hCorrCorrected[it][ia]->Add(hCorrSideMassProj[it][ia], -1);
            hCorrCorrected[it][ia]->Add(hCorrSideSideProj[it][ia]);

            hCorrCorrected[it][ia]->Scale(1.0/ncandidateswobg, "width");
        }
    }

}

void DoSimplifiedAnalysis()
{
    for (int it = 0; it < nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];
        int nreal = nRealTrigg[it];
        int ncandidates = hMassTrigg[it]->Integral(hMassTrigg[it]->GetXaxis()->FindBin(massMin), hMassTrigg[it]->GetXaxis()->FindBin(massMax));
        int ncandidateswobg = st[it];
        cout << "Trigg bin [" << tlow << " " << tupp << "] : \tN REAL TRIGGERS=" << nreal << "\tN CANDIDATES="  << ncandidates << "\tN CANDIDATES (NO BG)=" << ncandidateswobg << endl;
        for (int ia = 0; ia < nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            // True correlation
            hCorrRealProj[it][ia]->Scale(1.0/nreal, "width");

            // Non sideband corrected correlation
            //hCorrNonCorrected[it][ia]->Scale(1.0/fFitTrigg[it]->Integral(massMin, massMax), "width");
            hCorrNonCorrected[it][ia]->Scale(1.0/ncandidates, "width");

            // Sideband corrected correlation
            double scaling = (hCorrSideSideProj[it][ia]->GetEntries()+hCorrMassSideProj[it][ia]->GetEntries()+hCorrSideMassProj[it][ia]->GetEntries())/hCorrSideSideProj[it][ia]->GetEntries();
            hCorrSideSideProj[it][ia]->Scale(alpha[it][ia]*scaling);
            hCorrCorrected[it][ia]->Add(hCorrSideSideProj[it][ia], -1.);
            hCorrCorrected[it][ia]->Scale(1.0/ncandidateswobg, "width");
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

            if (tlow < aupp) continue;

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

            if (tlow < aupp) continue;

            // Assoc mass histrograms where trigger is from peak region
            fFitAssocPeak[it][ia] = new TF1(Form("fFitAssocPeak%d-%d", it, ia), FitFunction, 10, 500, 10);
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
        }
    }
}

void GetScaleFactorsVersion1()
{
    for (int it = 0; it < nTriggBins; it++) {
        st[it]  = hMassTrigg[it]->Integral(hMassTrigg[it]->GetXaxis()->FindBin(massMin), hMassTrigg[it]->GetXaxis()->FindBin(massMax)) - fBgTrigg[it]->Integral(massMin, massMax)/hMassTrigg[it]->GetBinWidth(0);
        //stobtrigg[it] = fPeakTrigg[it]->Integral(massMin, massMax)/fFitTrigg[it]->Integral(massMin, massMax);
        stobtrigg[it] = fPeakTrigg[it]->Integral(massMin, massMax)/fBgTrigg[it]->Integral(massMin, massMax);
        std::cout << "\n\tbin [ " << triggPt[it] << " " << triggPt[it+1] << " ] : \treal=" << nRealTrigg[it] << "\treconst=" << st[it] << "\trec=" << stfit[it] << "\trec/real=" << st[it]/nRealTrigg[it] << std::endl;
        //std::cout << "\tS/(S+B) : " << stobtrigg[it] << std::endl;
        std::cout << "\tS/B : " << stobtrigg[it] << std::endl;

        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;
            //alpha[it][ia]  = fBgAssocPeak[it][ia]->Integral(massMin, massMax)/(fBgAssocPeak[it][ia]->Integral(40, 80) + fBgAssocPeak[it][ia]->Integral(210, 280));
            //beeta[it][ia]  = fBgTrigg[it]->Integral(massMin, massMax)/(fBgTrigg[it]->Integral(40, 80) + fBgTrigg[it]->Integral(210, 280));
            //alpha[it][ia]  = fBgAssocPeak[it][ia]->Integral(massMin, massMax)/(fBgAssocPeak[it][ia]->Integral(50, 115) + fBgAssocPeak[it][ia]->Integral(160, 200));
            //beeta[it][ia]  = fBgTrigg[it]->Integral(massMin, massMax)/(fBgTrigg[it]->Integral(50, 115) + fBgTrigg[it]->Integral(160, 200));
            alpha[it][ia]  = fBgAssocPeak[it][ia]->Integral(massMin, massMax)/fBgAssocPeak[it][ia]->Integral(210, 350);
            beeta[it][ia]  = fBgTrigg[it]->Integral(massMin, massMax)/fBgTrigg[it]->Integral(210, 350);
            yamma[it][ia] = alpha[it][ia]*beeta[it][ia];
            std::cout << "\t\tbin [ "  << assocPt[ia] << " " << assocPt[ia+1] << " ] : "
                      << "\talpha=" << alpha[it][ia]
                      << "\tbeeta="  << beeta[it][ia]
                      << "\tgamma=" << yamma[it][ia] << std::endl;
            stobassoc[it][ia] = fPeakAssocPeak[it][ia]->Integral(massMin, massMax)/fFitAssocPeak[it][ia]->Integral(massMin, massMax);
            std::cout << "\t\tS/(S+B) : " << stobassoc[it][ia] << std::endl;
        }
    }
}

void DrawMassHistos(TString dataname)
{
    TLatex tl;
    tl.SetTextSize(0.035);
    tl.SetTextFont(42);
    for (int it = 0; it < nTriggBins; it++) {
        cMassTrigg[it] = new TCanvas(Form("cMassTrigg%d", it), Form("cMassTrigg%d", it), 600, 600);
        hMassTrigg[it]->SetTitle(Form("Trigger invariant mass, %.1f < p_{T,t} < %.1f GeV/c; M_{#gamma#gamma} (MeV/c^{2}); counts", triggPt[it], triggPt[it+1]));
        hMassTrigg[it]->GetXaxis()->SetRangeUser(0., 500.);
        hMassTrigg[it]->GetYaxis()->SetMaxDigits(3);

        fPeakColored[it] = (TF1*)fFitTrigg[it]->Clone("fPeakColored");
        fPeakColored[it]->SetFillColor(kRed-10);
        fPeakColored[it]->GetXaxis()->SetRangeUser(massMin, massMax);
        fPeakColored[it]->SetFillStyle(1001);
        fLeftSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fLeftSidebandColored[it]->SetFillColor(kGreen-10);
        fLeftSidebandColored[it]->GetXaxis()->SetRangeUser(40., 80.);
        fLeftSidebandColored[it]->SetFillStyle(1001);
        fRightSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fRightSidebandColored[it]->SetFillColor(kGreen-10);
        fRightSidebandColored[it]->GetXaxis()->SetRangeUser(210., 280.);
        fRightSidebandColored[it]->SetFillStyle(1001);
        fBgColored[it] = (TF1*)fBgTrigg[it]->Clone("fBgColored");
        fBgColored[it]->SetFillColor(kBlue-10);
        fBgColored[it]->GetXaxis()->SetRangeUser(massMin, massMax);
        fBgColored[it]->SetFillStyle(1001);

        cMassTrigg[it]->cd();
        hMassTrigg[it]->Draw("HIST");
        fPeakColored[it]->Draw("SAME");
        fLeftSidebandColored[it]->Draw("SAME");
        fRightSidebandColored[it]->Draw("SAME");
        fBgColored[it]->Draw("SAME");
        fFitTrigg[it]->Draw("SAME");
        fBgTrigg[it]->Draw("SAME");
        fPeakTrigg[it]->Draw("SAME");

        tl.DrawLatexNDC(.52, .7, "#pi^{0} mass peak");
        tl.DrawLatexNDC(.55, .65, Form("#bf{#bf{%.1f < p_{T,t} < %.1f GeV/c}}", triggPt[it], triggPt[it+1]));
        tl.DrawLatexNDC(.55, .6,  Form("S/B=%.3f", stobtrigg[it]));
        gPad->RedrawAxis();
        cMassTrigg[it]->SaveAs(Form("MassTrigg%d_%s.pdf", it, dataname.Data()));

        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            // Trigger from mass region
            cMassAssocPeak[it][ia] = new TCanvas(Form("cMassAssocPeak%d:%d",it,ia), Form("cMassAssocPeak%d:%d",it,ia), 600, 600);
            hMassAssocPeak[it][ia]->SetTitle("; M_{#gamma#gamma} (MeV/c^{2}); counts");
            //hMassAssocPeak[it][ia]->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; M_{#gamma#gamma}; counts", tlow, tupp, alow, aupp));
            hMassAssocPeak[it][ia]->GetXaxis()->SetRangeUser(0., 500.);
            hMassAssocPeak[it][ia]->GetYaxis()->SetMaxDigits(3);

            fPeakColoredAssocPeak[it][ia] = (TF1*)fFitAssocPeak[it][ia]->Clone("fPeakColoredAssocPeak");
            fPeakColoredAssocPeak[it][ia]->SetFillColor(kRed-10);
            fPeakColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(massMin, massMax);
            fPeakColoredAssocPeak[it][ia]->SetFillStyle(1001);
            fLeftSidebandColoredAssocPeak[it][ia] = (TF1*)fFitAssocPeak[it][ia]->Clone("fLeftSidebandColoredAssocPeak");
            fLeftSidebandColoredAssocPeak[it][ia]->SetFillColor(kGreen-10);
            fLeftSidebandColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(40., 80.);
            fLeftSidebandColoredAssocPeak[it][ia]->SetFillStyle(1001);
            fRightSidebandColoredAssocPeak[it][ia] = (TF1*)fFitAssocPeak[it][ia]->Clone("fLeftSidebandColoredAssocPeak");
            fRightSidebandColoredAssocPeak[it][ia]->SetFillColor(kGreen-10);
            fRightSidebandColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(210., 280.);
            fRightSidebandColoredAssocPeak[it][ia]->SetFillStyle(1001);
            fBgColoredAssocPeak[it][ia] = (TF1*)fBgAssocPeak[it][ia]->Clone("fBgColoredAssocPeak");
            fBgColoredAssocPeak[it][ia]->SetFillColor(kBlue-10);
            fBgColoredAssocPeak[it][ia]->GetXaxis()->SetRangeUser(massMin, massMax);
            fBgColoredAssocPeak[it][ia]->SetFillStyle(1001);

            cMassAssocPeak[it][ia]->cd();
            hMassAssocPeak[it][ia]->Draw("HIST");
            fPeakColoredAssocPeak[it][ia]->Draw("SAME");
            fLeftSidebandColoredAssocPeak[it][ia]->Draw("SAME");
            fRightSidebandColoredAssocPeak[it][ia]->Draw("SAME");
            fBgColoredAssocPeak[it][ia]->Draw("SAME");
            hMassAssocPeak[it][ia]->Draw("HIST SAME");
            fFitAssocPeak[it][ia]->Draw("SAME");
            fBgAssocPeak[it][ia]->Draw("SAME");
            fPeakAssocPeak[it][ia]->Draw("SAME");

            tl.DrawLatexNDC(.52, .4, "#pi^{0} mass peak");
            tl.DrawLatexNDC(.55, .35, Form("#bf{#bf{%.1f < p_{T,t} < %.1f GeV/c}}", tlow, tupp));
            tl.DrawLatexNDC(.55, .3, Form("%.1f < p_{T,a} < %.1f GeV/c", alow, aupp));
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
            fPeakColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(massMin, massMax);
            fPeakColoredAssocSide[it][ia]->SetFillStyle(1001);
            fLeftSidebandColoredAssocSide[it][ia] = (TF1*)fFitAssocSide[it][ia]->Clone("fLeftSidebandColoredAssocSide");
            fLeftSidebandColoredAssocSide[it][ia]->SetFillColor(kGreen-10);
            fLeftSidebandColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(40., 80.);
            fLeftSidebandColoredAssocSide[it][ia]->SetFillStyle(1001);
            fRightSidebandColoredAssocSide[it][ia] = (TF1*)fFitAssocSide[it][ia]->Clone("fLeftSidebandColoredAssocSide");
            fRightSidebandColoredAssocSide[it][ia]->SetFillColor(kGreen-10);
            fRightSidebandColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(210., 280.);
            fRightSidebandColoredAssocSide[it][ia]->SetFillStyle(1001);
            fBgColoredAssocSide[it][ia] = (TF1*)fBgAssocSide[it][ia]->Clone("fBgColoredAssocSide");
            fBgColoredAssocSide[it][ia]->SetFillColor(kBlue-10);
            fBgColoredAssocSide[it][ia]->GetXaxis()->SetRangeUser(massMin, massMax);
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
