#include "AnalyseCorrelations.h"

void AnalyseCorrelations()
{
    processDataSTAR();
    //processDataFoCal();
}

void processDataSTAR()
{
	TString fInName[ndata_star] = {
        "/home/heimarry/Simulations/focal/STAR-data/output_pAu.root",
        "/home/heimarry/Simulations/focal-pythia-sim/star-pp_WRONG-TRUE-CORR-COMP.root"
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
		GetScaleFactors();
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
        "/home/heimarry/Simulations/focal-pythia-sim/focal-all.root"
	};

	TString fOutName[ndata_focal] = {
		"analysis_FoCal_pp.root"
	};

    TString dataname[ndata_focal] = {
        "FoCal_pp"
    };

	for (int idata = 0; idata < ndata_focal; idata++) {
		fIn = TFile::Open(fInName[idata]);
		fOut = TFile::Open(fOutName[idata], "RECREATE");

		LoadInput();
        FitMassPeaks();
		GetScaleFactors();
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

	    for (int ia = 0; ia < nAssocBins; ia++) {
	        double alow = assocPt[ia];
	        double aupp = assocPt[ia+1];

	        if (tlow < aupp) continue;

	        hMassAssocPeak[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hMassAssocPeak[it][ia]->Rebin();
	        hMassAssocSide[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); //hMassAssocSide[it][ia]->Rebin();

	        hCorrReal[it][ia]     = (TH2D*)fIn->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrReal[it][ia]->Rebin2D(4);
	        hCorrMassMass[it][ia] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrMassMass[it][ia]->Rebin2D(4);
	        hCorrMassSide[it][ia] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrMassSide[it][ia]->Rebin2D(4);
	        hCorrSideMass[it][ia] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrSideMass[it][ia]->Rebin2D(4);
	        hCorrSideSide[it][ia] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrSideSide[it][ia]->Rebin2D(4);
            hCorrMeas[it][ia]     = (TH2D*)hCorrMassMass[it][ia]->Clone(Form("hCorrMeas[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

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

            if (tlow < aupp) continue;

            hCorrRealProj[it][ia] = hCorrReal[it][ia]->ProjectionX();
            hCorrRealProj[it][ia]->Scale(1.0/nRealTrigg[it], "width");

            hCorrMassMassProj[it][ia] = hCorrMassMass[it][ia]->ProjectionX();

            hCorrMassSideProj[it][ia] = hCorrMassSide[it][ia]->ProjectionX();
            hCorrMassSideProj[it][ia]->Scale(alpha[it][ia]);

            hCorrSideMassProj[it][ia] = hCorrSideMass[it][ia]->ProjectionX();
            hCorrSideMassProj[it][ia]->Scale(beta[it]);

            hCorrSideSideProj[it][ia] = hCorrSideSide[it][ia]->ProjectionX();
            hCorrSideSideProj[it][ia]->Scale(alpha[it][ia]*beta[it]);

            hCorr[it][ia] = (TH1D*)hCorrMassMassProj[it][ia]->Clone(Form("hCorrFinal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorr[it][ia]->Add(hCorrMassSideProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideMassProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideSideProj[it][ia]);
            hCorr[it][ia]->Scale(1.0/st[it], "width");

            hCorrMeasProj[it][ia] = hCorrMeas[it][ia]->ProjectionX();
            hCorrMeasProj[it][ia]->Scale(1.0/st[it], "width");
        }
    }

}

void FitMassPeaks()
{
    for (int it = 0; it < nTriggBins; it++) {
        fFitTrigg[it] = new TF1(Form("fFitTrigg%d", it), FitFunction, 20, 300, 6);
        fFitTrigg[it]->SetParameters(1., 1., 1., 135., 10., 15.);
        fFitTrigg[it]->SetParLimits(3, 132., 137.);
        fFitTrigg[it]->SetNpx(1000);

        fPeakTrigg[it] = new TF1(Form("fPeakTrigg%d", it), FitPeak, 20, 300, 5);
        fPeakTrigg[it]->SetLineColor(kBlue);
        fPeakTrigg[it]->SetNpx(1000);

        fBgTrigg[it] = new TF1(Form("fBgTrigg%d", it), FitBackground, 20, 300, 4);
        fBgTrigg[it]->SetLineColor(kBlack);
        fBgTrigg[it]->SetLineStyle(kDashed);
        fBgTrigg[it]->SetNpx(1000);

        hMassTrigg[it]->Fit(Form("fFitTrigg%d", it), "QNL0R+");
        fFitTrigg[it]->GetParameters(parTrigg);
        fBgTrigg[it]->SetParameters(parTrigg);
        fPeakTrigg[it]->SetParameters(&parTrigg[3]);

        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            fFitAssoc[it][ia] = new TF1(Form("fFitAssoc%d-%d", it, ia), FitFunction, 10, 300, 6);
            fFitAssoc[it][ia]->SetParameters(1., 1., 1., 135., 30., 50.);
            fFitAssoc[it][ia]->SetParLimits(3, 132., 137.);
            fFitAssoc[it][ia]->SetNpx(1000);

            fPeakAssoc[it][ia] = new TF1(Form("fPeakAssoc%d-%d", it, ia), FitPeak, 10, 300, 5);
            fPeakAssoc[it][ia]->SetLineColor(kBlue);
            fPeakAssoc[it][ia]->SetNpx(1000);

            fBgAssoc[it][ia]= new TF1(Form("fBgAssoc%d-%d", it, ia), FitBackground, 10, 300, 3);
            fBgAssoc[it][ia]->SetLineColor(kBlack);
            fBgAssoc[it][ia]->SetLineStyle(kDashed);
            fBgAssoc[it][ia]->SetNpx(1000);

            hMassAssocPeak[it][ia]->Fit(Form("fFitAssoc%d-%d", it, ia), "QNL0R+");
            fFitAssoc[it][ia]->GetParameters(parAssoc);
            fBgAssoc[it][ia]->SetParameters(parAssoc);
            fPeakAssoc[it][ia]->SetParameters(&parAssoc[3]);
        }
    }
}

void GetScaleFactors()
{
    for (int it = 0; it < nTriggBins; it++) {
        beta[it] = fBgTrigg[it]->Integral(110, 160)/(fBgTrigg[it]->Integral(40, 80) + fBgTrigg[it]->Integral(210, 280));
        st[it] = fPeakTrigg[it]->Integral(110, 160);
        std::cout << "\n\tbin [ " << triggPt[it] << " " << triggPt[it+1] << " ] : \treal=" << nRealTrigg[it] << "\treconst=" << st[it] << "\trec/real=" << st[it]/nRealTrigg[it] << std::endl;
        std::cout << "\tbin [ " << triggPt[it] << " " << triggPt[it+1] << " ] : \tbeta=" << beta[it] << std::endl;
        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;
            alpha[it][ia] = fBgAssoc[it][ia]->Integral(110, 160)/(fBgAssoc[it][ia]->Integral(40, 80) + fBgAssoc[it][ia]->Integral(210, 280));
            std::cout << "\t\tbin [ " << assocPt[ia] << " " << assocPt[ia+1] << " ]" << " : \talpha=" << alpha[it][ia] << std::endl;
        }
    }
}

void DrawMassHistos(TString dataname)
{
    for (int it = 0; it < nTriggBins; it++) {
        cMassTrigg[it] = new TCanvas(Form("cMassTrigg%d", it), Form("cMassTrigg%d", it), 600, 600);
        hMassTrigg[it]->SetTitle(Form("Trigger invariant mass, %.1f < p_{T,t} < %.1f GeV/c; M_{#gamma#gamma} (MeV/c^{2}); counts", triggPt[it], triggPt[it+1]));
        hMassTrigg[it]->GetXaxis()->SetRangeUser(0., 300.);
        hMassTrigg[it]->GetYaxis()->SetMaxDigits(3);

        fPeakColored[it] = (TF1*)fFitTrigg[it]->Clone("fPeakColored");
        fPeakColored[it]->SetFillColor(kBlue-10);
        fPeakColored[it]->GetXaxis()->SetRangeUser(110., 160.);
        fPeakColored[it]->SetFillStyle(1001);
        fLeftSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fLeftSidebandColored[it]->SetFillColor(kRed-10);
        fLeftSidebandColored[it]->GetXaxis()->SetRangeUser(40., 80.);
        fLeftSidebandColored[it]->SetFillStyle(1001);
        fRightSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fRightSidebandColored[it]->SetFillColor(kRed-10);
        fRightSidebandColored[it]->GetXaxis()->SetRangeUser(210., 280.);
        fRightSidebandColored[it]->SetFillStyle(1001);

        cMassTrigg[it]->cd();
        hMassTrigg[it]->Draw("HIST");
        fPeakColored[it]->Draw("SAME");
        fLeftSidebandColored[it]->Draw("SAME");
        fRightSidebandColored[it]->Draw("SAME");
        fFitTrigg[it]->Draw("SAME");
        fBgTrigg[it]->Draw("SAME");
        fPeakTrigg[it]->Draw("SAME");
        gPad->RedrawAxis();
        cMassTrigg[it]->SaveAs(Form("MassTrigg%d_%s.pdf", it, dataname.Data()));

        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            cMassAssocPeak[it][ia] = new TCanvas(Form("cMassAssocPeak%d:%d",it,ia), Form("cMassAssocPeak%d:%d",it,ia), 600, 600);
            hMassAssocPeak[it][ia]->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; M_{#gamma#gamma}; counts", tlow, tupp, alow, aupp));
            hMassAssocPeak[it][ia]->GetXaxis()->SetRangeUser(0., 300.);
            hMassAssocPeak[it][ia]->GetYaxis()->SetMaxDigits(6);

            fPeakColoredAssoc[it][ia] = (TF1*)fFitAssoc[it][ia]->Clone("fPeakColoredAssoc");
            fPeakColoredAssoc[it][ia]->SetFillColor(kBlue-10);
            fPeakColoredAssoc[it][ia]->GetXaxis()->SetRangeUser(110., 160.);
            fPeakColoredAssoc[it][ia]->SetFillStyle(1001);
            fLeftSidebandColoredAssoc[it][ia] = (TF1*)fFitAssoc[it][ia]->Clone("fLeftSidebandColoredAssoc");
            fLeftSidebandColoredAssoc[it][ia]->SetFillColor(kRed-10);
            fLeftSidebandColoredAssoc[it][ia]->GetXaxis()->SetRangeUser(40., 80.);
            fLeftSidebandColoredAssoc[it][ia]->SetFillStyle(1001);
            fRightSidebandColoredAssoc[it][ia] = (TF1*)fFitAssoc[it][ia]->Clone("fLeftSidebandColoredAssoc");
            fRightSidebandColoredAssoc[it][ia]->SetFillColor(kRed-10);
            fRightSidebandColoredAssoc[it][ia]->GetXaxis()->SetRangeUser(210., 280.);
            fRightSidebandColoredAssoc[it][ia]->SetFillStyle(1001);

            cMassAssocPeak[it][ia]->cd();
            hMassAssocPeak[it][ia]->Draw("HIST");
            fPeakColoredAssoc[it][ia]->Draw("SAME");
            fLeftSidebandColoredAssoc[it][ia]->Draw("SAME");
            fRightSidebandColoredAssoc[it][ia]->Draw("SAME");
            fFitAssoc[it][ia]->Draw("SAME");
            fBgAssoc[it][ia]->Draw("SAME");
            fPeakAssoc[it][ia]->Draw("SAME");
            gPad->RedrawAxis();
            cMassAssocPeak[it][ia]->SaveAs(Form("MassAssocPeak%d-%d_%s.pdf",it,ia,dataname.Data()));
        }
    }
}

double FitPeak(double *x, double *p)
{
    double c1 = p[0];
    //double c2 = p[1];
    double mu = p[1];
    double sigma1 = p[2];
    //double sigma2 = p[4];

    return c1*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma1*sigma1));// + c2*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma2*sigma2));
}

//double FitPeak(double *x, double *p) {
//    double p1 = p[0];
//    double p2 = p[1];
//    double p3 = p[2];
//    return (0.5*p[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
//}

double FitBackground(double *x, double *p)
{
    double a = p[0];
    double b = p[1];
    double c = p[2];

    return a*x[0]*x[0] + b*x[0] + c;
}

double FitFunction(double *x, double *p)
{
    return FitBackground(x, p) + FitPeak(x, &p[3]);
}
