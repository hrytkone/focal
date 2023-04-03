#include "include/Filipad.h"
#include "include/rootcommon.h"

void LoadData();
void ConfigHistos();
void DrawFiliPad();

const int nTriggBins = 2;
const int nAssocBins = 2;

const double triggPt[nTriggBins+1] = {4.0, 8.0, 16.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};
//const double triggPt[nTriggBins+1] = {1.0, 2.0, 2.5, 3.0};
//const double assocPt[nAssocBins+1] = {0.5, 1.0, 1.5, 2.0, 2.5};
const double mSize = 1.;

const double ymin[nTriggBins][nAssocBins] = {{0.008, 0.0008}, {0.0001, 0.0001}};
const double ymax[nTriggBins][nAssocBins] = {{2, 2}, {2, 2}};

const int nset = 1;

TString infiles[nset] = {
    //"analysis_FoCal_pp.root"
    //"analysis_FoCal_pp_pythia.root"
    //"analysis_FoCal_pPb_pythia.root"
    "analysis_FoCal_pp_geant.root"
};

TString legHeader[nset] = {
    "p-p #sqrt{s} = 14 TeV"
    //"pPb #sqrt{s} = 5.02 TeV"
};

TFile *fin[nset];

Filipad *fpad[nset][nTriggBins][nAssocBins];
Filipad *fpadPtComp;

TH1D *hCorrReal[nset][nTriggBins][nAssocBins];
TH1D *hCorrCorrected[nset][nTriggBins][nAssocBins];
TH1D *hCorrNonCorrected[nset][nTriggBins][nAssocBins];
TH1D *hRatioCorrected[nset][nTriggBins][nAssocBins];
TH1D *hRatioNonCorrected[nset][nTriggBins][nAssocBins];

TLegend *leg[nset][nTriggBins][nAssocBins];

void ComparisonPlot()
{
    gStyle->SetOptStat(0);

	LoadData();
    ConfigHistos();
	DrawFiliPad();
}

void LoadData()
{
	for (int iset = 0; iset < nset; iset++)
		fin[iset] = TFile::Open(infiles[iset]);

	for (int iset = 0; iset < nset; iset++) {
		for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
			for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
    	        double aupp = assocPt[iassoc+1];

    	        if (tlow < aupp) continue;

                hCorrReal[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hCorrCorrected[iset][itrigg][iassoc] = (TH1D*)fin[iset]->Get(Form("hCorrCorrected[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hCorrNonCorrected[iset][itrigg][iassoc] = (TH1D*)fin[iset]->Get(Form("hCorrNonCorrected[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

                cout << "BIN WIDTH : " << hCorrReal[iset][itrigg][iassoc]->GetBinWidth(0) << endl;

                //// Rebin if needed
                //hCorrReal[iset][itrigg][iassoc]->Rebin(8);
                //hCorrCorrected[iset][itrigg][iassoc]->Rebin(8);
                //hCorrNonCorrected[iset][itrigg][iassoc]->Rebin(8);

                // Calculate ratios
                hRatioCorrected[iset][itrigg][iassoc] = (TH1D*)hCorrCorrected[iset][itrigg][iassoc]->Clone(Form("hRatioCorrected[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioCorrected[iset][itrigg][iassoc]->Divide(hCorrReal[iset][itrigg][iassoc]);
                hRatioNonCorrected[iset][itrigg][iassoc] = (TH1D*)hCorrNonCorrected[iset][itrigg][iassoc]->Clone(Form("hRatioNonCorrected[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioNonCorrected[iset][itrigg][iassoc]->Divide(hCorrReal[iset][itrigg][iassoc]);
			}
		}
	}
}

void ConfigHistos()
{
    for (int iset = 0; iset < nset; iset++) {
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;

                hCorrReal[iset][itrigg][iassoc]->SetLineColor(kBlack);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerStyle(kOpenSquare);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerSize(mSize);
                hCorrCorrected[iset][itrigg][iassoc]->SetLineColor(kRed);
                hCorrCorrected[iset][itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrCorrected[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hCorrCorrected[iset][itrigg][iassoc]->SetMarkerSize(mSize);
                hCorrNonCorrected[iset][itrigg][iassoc]->SetLineColor(kBlue);
                hCorrNonCorrected[iset][itrigg][iassoc]->SetMarkerColor(kBlue);
                hCorrNonCorrected[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hCorrNonCorrected[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hRatioCorrected[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(0., 5.);
                hRatioCorrected[iset][itrigg][iassoc]->SetLineColor(kRed);
                hRatioCorrected[iset][itrigg][iassoc]->SetMarkerColor(kRed);
                hRatioCorrected[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hRatioCorrected[iset][itrigg][iassoc]->SetMarkerSize(mSize);
                hRatioNonCorrected[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(0., 5.);
                hRatioNonCorrected[iset][itrigg][iassoc]->SetLineColor(kBlue);
                hRatioNonCorrected[iset][itrigg][iassoc]->SetMarkerColor(kBlue);
                hRatioNonCorrected[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hRatioNonCorrected[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                //leg[iset][itrigg][iassoc] = new TLegend(0.58, 0.45, 0.85, 0.75);
                leg[iset][itrigg][iassoc] = new TLegend(0.58, 0.4, 0.85, 0.7);
                leg[iset][itrigg][iassoc]->SetFillStyle(0); leg[iset][itrigg][iassoc]->SetBorderSize(0); leg[iset][itrigg][iassoc]->SetTextSize(0.05);
                leg[iset][itrigg][iassoc]->SetHeader(Form("#splitline{%s}{[%0.1f,%0.1f][%0.1f,%0.1f]}", legHeader[iset].Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
                leg[iset][itrigg][iassoc]->AddEntry(hCorrReal[0][0][0], "MC truth", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrNonCorrected[0][0][0], "Before SB corr", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrCorrected[0][0][0], "After SB corr", "pe");
            }
        }
    }
}

void DrawFiliPad()
{
    TText *t = new TText(.58,0.79,"PYTHIA8 simulation");
    //TText *t = new TText(.58,0.79," PYTHIA6+GEANT3");
    t->SetNDC();

    TLatex tl;
    tl.SetTextSize(0.05);
    //tl.SetTextFont(42);

    int padID = 0;
    for (int iset = 0; iset < nset; iset++) {
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;

//                fpad[iset][itrigg][iassoc] = new Filipad(padID, 1.1, 0.5, 100, 100, 0.7, 1);
                fpad[iset][itrigg][iassoc] = new Filipad(padID, 1.1, 0.35, 100, 100, 0.8, 1);
                fpad[iset][itrigg][iassoc]->Draw();
                padID++;

                // Upper pad
                int minBin = hCorrCorrected[iset][itrigg][iassoc]->GetMinimumBin();
                int maxBin = hCorrNonCorrected[iset][itrigg][iassoc]->GetMaximumBin();
                double rangeMin = hCorrCorrected[iset][itrigg][iassoc]->GetBinContent(minBin) - 0.1*hCorrCorrected[iset][itrigg][iassoc]->GetBinContent(minBin);
                double rangeMax = hCorrNonCorrected[iset][itrigg][iassoc]->GetBinContent(maxBin) + 10.*hCorrNonCorrected[iset][itrigg][iassoc]->GetBinContent(maxBin);
                if (rangeMin<=0) rangeMin = 0.0008;
                rangeMin = ymin[itrigg][iassoc];
                rangeMax = ymax[itrigg][iassoc];
                cout << rangeMin << " " << rangeMax << endl;
                TPad *p = fpad[iset][itrigg][iassoc]->GetPad(1);
                p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
                hset(*hCorrNonCorrected[iset][itrigg][iassoc], "#Delta#phi", "1/N_{trigg}dN/d#Delta#phi", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
//                hset(*hCorrReal[iset][itrigg][iassoc], "#Delta#phi", "arb. norm.", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
                hCorrNonCorrected[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(rangeMin, rangeMax);
//                hCorrReal[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(TMath::Pi()/2., 3.*TMath::Pi()/2.);
                hCorrNonCorrected[iset][itrigg][iassoc]->Draw("P");
                hCorrCorrected[iset][itrigg][iassoc]->Draw("SAME P");
                hCorrReal[iset][itrigg][iassoc]->Draw("SAME P");
                leg[iset][itrigg][iassoc]->Draw("SAME");
                t->Draw("SAME");
                tl.DrawLatexNDC(0.58, 0.74, "Asymmetry < 0.8");
                //tl.DrawLatexNDC(0.58, 0.74, "No asymmetry cut");

                // Lower pad
                p = fpad[iset][itrigg][iassoc]->GetPad(2);
                p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
                hset( *hRatioNonCorrected[iset][itrigg][iassoc], "#Delta#phi", "Ratio to MC truth",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
                hRatioNonCorrected[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.2, 1.8);
                //hRatioNonCorrected[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(TMath::Pi()/2., 3.*TMath::Pi()/2.);
                hRatioNonCorrected[iset][itrigg][iassoc]->Draw("P");
                hRatioCorrected[iset][itrigg][iassoc]->Draw("P SAME");
                //hRatioSBMeas[iset][itrigg][iassoc]->Draw("P SAME");
            }
        }
    }
}
