#include "AnalyseCorrelations.h"

void AnalyseCorrelations()
{
    processDataSTAR();
}

void processDataSTAR()
{
	TString fInName[ndata] = {
        "/mnt/STAR-sim/output_pAu.root",
        "/mnt/STAR-sim/output_pp.root"
	};

	TString fOutName[ndata] = {
        "analysis_STAR_pAu.root",
		"analysis_STAR_pp.root"
	};


	for (int idata = 0; idata < ndata; idata++) {
		fIn = TFile::Open(fInName[idata]);
		fOut = TFile::Open(fOutName[idata], "RECREATE");

		LoadInput();
		GetNormalisations();
		DoAnalysis();
        PlotMassHistos();
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
	    //hMassTrigg[it]->Scale(1./nEvent);
	    //hMassTrigg[it]->Rebin();

	    for (int ia = 0; ia < nAssocBins; ia++) {
	        double alow = assocPt[ia];
	        double aupp = assocPt[ia+1];

	        if (tlow < aupp) continue;

	        hMassAssocPeak[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
	        //hMassAssocPeak[it][ia]->Scale(1./nEvent);

	        hMassAssocSide[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
	        //hMassAssocSide[it][ia]->Scale(1./nEvent);

	        hCorrReal[it][ia] = (TH2D*)fIn->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrReal[it][ia]->Rebin2D(4);
	        hCorrMassMass[it][ia] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrMassMass[it][ia]->Rebin2D(4);
	        hCorrMassSide[it][ia] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrMassSide[it][ia]->Rebin2D(4);
	        hCorrSideMass[it][ia] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrSideMass[it][ia]->Rebin2D(4);
	        hCorrSideSide[it][ia] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrSideSide[it][ia]->Rebin2D(4);
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

            double nPairReal = hCorrReal[it][ia]->GetEntries();
            double nPairMassMass = hCorrMassMass[it][ia]->GetEntries();
            double nPairMassSide = hCorrMassSide[it][ia]->GetEntries();
            double nPairSideMass = hCorrSideMass[it][ia]->GetEntries();
            double nPairSideSide = hCorrSideSide[it][ia]->GetEntries();

            hCorrRealProj[it][ia] = hCorrReal[it][ia]->ProjectionX();
            //hCorrRealProj[it][ia]->Scale(1.0/nPairMassMass);

            hCorrMassMassProj[it][ia] = hCorrMassMass[it][ia]->ProjectionX();
            //hCorrMassMassProj[it][ia]->Scale(1.0/nPairMassMass);
            //hCorrMassMassProj[it][ia]->Scale(1/(alphat[it]*alphaa[it][ia]));

            hCorrMassSideProj[it][ia] = hCorrMassSide[it][ia]->ProjectionX();
            //hCorrMassSideProj[it][ia]->Scale(1.0/nPairMassSide);
            //hCorrMassSideProj[it][ia]->Scale(betaa[it][ia]/(alphat[it]*alphaa[it][ia]));
            //hCorrMassSideProj[it][ia]->Scale(B[it][ia]);

            hCorrSideMassProj[it][ia] = hCorrSideMass[it][ia]->ProjectionX();
            //hCorrSideMassProj[it][ia]->Scale(1.0/nPairSideMass);
            //hCorrSideMassProj[it][ia]->Scale(betat[it]/(alphat[it]*alphaa[it][ia]));
            //hCorrSideMassProj[it][ia]->Scale(A[it][ia]);

            hCorrSideSideProj[it][ia] = hCorrSideSide[it][ia]->ProjectionX();
            //hCorrSideSideProj[it][ia]->Scale(1.0/nPairSideSide);
            //hCorrSideSideProj[it][ia]->Scale((betat[it]*betaa[it][ia])/(alphat[it]*alphaa[it][ia]));
            //hCorrSideSideProj[it][ia]->Scale((A[it][ia]*B[it][ia]));

            hCorr[it][ia] = (TH1D*)hCorrMassMassProj[it][ia]->Clone(Form("hCorr%d:%d",it,ia));
            hCorr[it][ia]->Add(hCorrMassSideProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideMassProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideSideProj[it][ia]);
            //hCorr[it][ia]->Scale(1./(brGammaCh*brGammaCh)); // Correct for branching ratio
            //hCorr[it][ia]->Scale(1./(0.95*0.95)); // Correct for missing photon pair at limit of acceptance
        }
    }

}

void GetNormalisations()
{
	// Mass peak fits
	double parTrigg[8];
	double parAssoc[9];
    for (int it = 0; it < nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];

        //fFitTrigg[it] = new TF1(Form("fFitTrigg[%4.1f,%4.1f]",tlow,tupp), FitFunction, 20, 300, 8);
        fFitTrigg[it] = new TF1(Form("fFitTrigg[%4.1f,%4.1f]",tlow,tupp), "[0]*TMath::Exp(-(x[0]-[2])*(x[0]-[2])/(2*[3]*[3])) + [1]*TMath::Exp(-(x[0]-[2])*(x[0]-[2])/(2*[4]*[4]))", 20, 300);
        fFitTrigg[it]->SetParameters(-1., 1., 1., 100., 100., 135., 10., 15.);
        fFitTrigg[it]->SetParLimits(5, 132, 138);
        fFitTrigg[it]->SetParLimits(6, 1., 15.);
        fFitTrigg[it]->SetParLimits(7, 2., 20.);
        fFitTrigg[it]->SetNpx(1000);

        fPeakTrigg[it] = new TF1(Form("fPeakTrigg[%4.1f,%4.1f]",tlow,tupp), FitPeak, 20, 300, 5);
        fPeakTrigg[it]->SetNpx(1000);

        fBgTrigg[it] = new TF1(Form("fBgTrigg[%4.1f,%4.1f]",tlow,tupp), FitBackground, 20, 300, 4);
        fBgTrigg[it]->SetNpx(1000);

        hMassTrigg[it]->Fit(Form("fFitTrigg[%4.1f,%4.1f]",tlow,tupp), "Q0R+");
        fFitTrigg[it]->GetParameters(parTrigg);
        fBgTrigg[it]->SetParameters(parTrigg);
        fPeakTrigg[it]->SetParameters(&parTrigg[3]);

        double sumt = fFitTrigg[it]->Integral(110, 160);
        double st = fPeakTrigg[it]->Integral(110, 160);
        double bt = fBgTrigg[it]->Integral(110, 160);
        double bsidet = fBgTrigg[it]->Integral(40, 80) + fBgTrigg[it]->Integral(210, 280);
        alphat[it] = st/sumt;
        betat[it] = bt/sumt;

        for (int ia = 0; ia < nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            fFitAssoc[it][ia] = new TF1(Form("fFitAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), FitFunction, 10, 300, 8);
            fFitAssoc[it][ia]->SetParameters(0., 0., 0., 0., 0., 135., 30., 50.);
            fFitAssoc[it][ia]->SetParLimits(5, 132, 138);
            fFitAssoc[it][ia]->SetParLimits(6, 2., 15.);
            fFitAssoc[it][ia]->SetParLimits(7, 4., 20.);

            fFitAssoc[it][ia]->SetNpx(1000);

            fPeakAssoc[it][ia] = new TF1(Form("fPeakAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), FitPeak, 10, 300, 5);
            fPeakAssoc[it][ia]->SetNpx(1000);

            fBgAssoc[it][ia]= new TF1(Form("fBgAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), FitBackground, 10, 300, 3);
            fBgAssoc[it][ia]->SetNpx(1000);

            hMassAssocPeak[it][ia]->Fit(Form("fFitAssoc[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp), "Q0R+");
            fFitAssoc[it][ia]->GetParameters(parAssoc);
            fBgAssoc[it][ia]->SetParameters(parAssoc);
            fPeakAssoc[it][ia]->SetParameters(&parAssoc[3]);

            double suma = fFitAssoc[it][ia]->Integral(110, 160);
            double sa = fPeakAssoc[it][ia]->Integral(110, 160);
            double ba = fBgAssoc[it][ia]->Integral(110, 160);
            double bsidea = fBgAssoc[it][ia]->Integral(40, 80) + fBgAssoc[it][ia]->Integral(210, 280);
            alphaa[it][ia] = sa/suma;
            betaa[it][ia] = ba/suma;
            A[it][ia] = bt/bsidet;
            B[it][ia] = ba/bsidea;

            //cout << "Normalization factors : A = " << A[it][ia] << "\tB = " << B[it][ia] << endl;
        }
    }
}

void PlotMassHistos()
{
    for (int i=0; i<ndata; i++) {
        cMassesTrigg[i] = new TCanvas(Form("cMassesTrigg%d",i), "", 800, 400);
        cMassesTrigg[i]->Divide(2,1);
        for (int it=0; it<nTriggBins; it++) {
            cMassesTrigg[i]->cd(it+1);
            hMassTrigg[it]->Draw();
        }
        cMassesTrigg[i]->SaveAs(Form("masstrigg%d.pdf",i));
    }
}

double FitPeak(double *x, double *p)
{
    double c1 = p[0];
    double c2 = p[1];
    double mu = p[2];
    double sigma1 = p[3];
    double sigma2 = p[4];

    return c1*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma1*sigma1)) + c2*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma2*sigma2));
}

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
