#include "include/Filipad.h"
#include "include/rootcommon.h"

void LoadData();
void ConfigHistos();
void DrawFiliPad();

const bool useLeading = 0;

const int nTriggBins = 2;
const int nAssocBins = 2;

const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};
//const double triggPt[nTriggBins+1] = {8.0, 20.0};
//const double assocPt[nAssocBins+1] = {3.0, 4.0, 8.0};
//const double triggPt[nTriggBins+1] = {4.0, 8.0, 10.0, 15.0, 20.0};
//const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0, 8.0, 10.0, 15.0};

// LEADING PARTICLE
//const double triggPt[nTriggBins+1] = {4.0, 10000.0};;
//const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};;

const double mSize = 1.;

const int nset = 1;
//TString infiles[nset] = {
//    "analysis_STAR_pp.root"//,
//	//"analysis_STAR_pAu.root"
//};

TString infiles[nset] = {
    //"analysis_FoCal_pp_fullsim_no-mixed.root"
    "analysis_FoCal_pp.root"
    //"analysis_FoCal_pp_check.root"
    //"analysis_FoCal_pp_fullsim_mixed.root",
    //"analysis_FoCal_pp_test-pythia.root"
};

TString legHeader[nset] = {
    "p-p #sqrt{s} = 14 TeV"
    //"p-p #sqrt{s} = 14 TeV"//,
    //"p-Au"
};

TFile *fin[nset];
Filipad *fpad[nset][nTriggBins][nAssocBins];
Filipad *fpadPtComp;
TH1D *hCorrReal[nset][nTriggBins][nAssocBins];
TH1D *hCorrFinal[nset][nTriggBins][nAssocBins];
TH1D *hCorrMeas[nset][nTriggBins][nAssocBins];
TH1D *hRatioFinal[nset][nTriggBins][nAssocBins];
TH1D *hRatioMeas[nset][nTriggBins][nAssocBins];
TH1D *hRatioSBMeas[nset][nTriggBins][nAssocBins];
TH1D *hRatioPt;
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

    	        if (!useLeading && tlow < aupp) continue;

                hCorrReal[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));  hCorrReal[iset][itrigg][iassoc]->Rebin(4);
                hCorrFinal[iset][itrigg][iassoc] = (TH1D*)fin[iset]->Get(Form("hCorrFinal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));   hCorrFinal[iset][itrigg][iassoc]->Rebin(4);
                hCorrMeas[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrMeas[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp)); hCorrMeas[iset][itrigg][iassoc]->Rebin(4);

                // Calculate ratios
                hRatioFinal[iset][itrigg][iassoc] = (TH1D*)hCorrReal[iset][itrigg][iassoc]->Clone(Form("hRatioFinal[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioFinal[iset][itrigg][iassoc]->Divide(hCorrFinal[iset][itrigg][iassoc]);
                hRatioMeas[iset][itrigg][iassoc] = (TH1D*)hCorrReal[iset][itrigg][iassoc]->Clone(Form("hRatioMeas[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioMeas[iset][itrigg][iassoc]->Divide(hCorrMeas[iset][itrigg][iassoc]);
                hRatioSBMeas[iset][itrigg][iassoc] = (TH1D*)hCorrFinal[iset][itrigg][iassoc]->Clone(Form("hRatioSBMeas[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioSBMeas[iset][itrigg][iassoc]->Divide(hCorrMeas[iset][itrigg][iassoc]);
			}
		}
	}
    hRatioPt = (TH1D*)hCorrReal[0][1][1]->Clone("hRatioPt");
    hRatioPt->Divide(hCorrReal[0][0][0]);
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

                if (!useLeading && tlow < aupp) continue;

                hCorrReal[iset][itrigg][iassoc]->SetLineColor(kBlack);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerSize(mSize);
                hCorrFinal[iset][itrigg][iassoc]->SetLineColor(kRed);
                hCorrFinal[iset][itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrFinal[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hCorrFinal[iset][itrigg][iassoc]->SetMarkerSize(mSize);
                if (iset==0) {
                    hCorrMeas[iset][itrigg][iassoc]->SetLineColor(kBlue);
                    hCorrMeas[iset][itrigg][iassoc]->SetMarkerColor(kBlue);
                } else {
                    hCorrMeas[iset][itrigg][iassoc]->SetLineColor(kRed);
                    hCorrMeas[iset][itrigg][iassoc]->SetMarkerColor(kRed);
                }

                hCorrMeas[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hCorrMeas[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hRatioFinal[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(0., 5.);
                hRatioFinal[iset][itrigg][iassoc]->SetLineColor(kRed-9);
                hRatioFinal[iset][itrigg][iassoc]->SetMarkerColor(kRed-9);
                hRatioFinal[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hRatioFinal[iset][itrigg][iassoc]->SetMarkerSize(mSize);
                hRatioMeas[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(0., 5.);
                if (iset==0) {
                    hRatioMeas[iset][itrigg][iassoc]->SetLineColor(kBlue-9);
                    hRatioMeas[iset][itrigg][iassoc]->SetMarkerColor(kBlue-9);
                } else {
                    hRatioMeas[iset][itrigg][iassoc]->SetLineColor(kRed-9);
                    hRatioMeas[iset][itrigg][iassoc]->SetMarkerColor(kRed-9);
                }

                hRatioMeas[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hRatioMeas[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hRatioSBMeas[iset][itrigg][iassoc]->SetLineColor(kViolet);
                hRatioSBMeas[iset][itrigg][iassoc]->SetMarkerColor(kViolet);
                hRatioSBMeas[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hRatioSBMeas[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                leg[iset][itrigg][iassoc] = new TLegend(0.58, 0.45, 0.85, 0.75);
                leg[iset][itrigg][iassoc]->SetFillStyle(0); leg[iset][itrigg][iassoc]->SetBorderSize(0); leg[iset][itrigg][iassoc]->SetTextSize(0.05);
                leg[iset][itrigg][iassoc]->SetHeader(Form("#splitline{%s}{[%0.1f,%0.1f][%0.1f,%0.1f]}", legHeader[iset].Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
                leg[iset][itrigg][iassoc]->AddEntry(hCorrReal[0][0][0], "MC truth", "pe");
                //leg[iset][itrigg][iassoc]->AddEntry(hCorrMeas[0][0][0], "GEANT3 sim", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrMeas[0][0][0], "Before SB corr", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrFinal[0][0][0], "After SB corr", "pe");
            }
        }
    }
}

void DrawFiliPad()
{
    TText *t = new TText(.58,0.79,"PYTHIA8 simulation");
    //TText *t = new TText(.58,0.79," FOCAL simulation");
    t->SetNDC();
    int padID = 0;
    for (int iset = 0; iset < nset; iset++) {
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (!useLeading && tlow < aupp) continue;

//                fpad[iset][itrigg][iassoc] = new Filipad(padID, 1.1, 0.5, 100, 100, 0.7, 1);
                fpad[iset][itrigg][iassoc] = new Filipad(padID, 1.1, 0.35, 100, 100, 0.8, 1);
                fpad[iset][itrigg][iassoc]->Draw();
                padID++;

                // Upper pad
                int minBin = hCorrFinal[iset][itrigg][iassoc]->GetMinimumBin();
                int maxBin = hCorrMeas[iset][itrigg][iassoc]->GetMaximumBin();
                double rangeMin = hCorrFinal[iset][itrigg][iassoc]->GetBinContent(minBin) - 0.1*hCorrFinal[iset][itrigg][iassoc]->GetBinContent(minBin);
                double rangeMax = hCorrMeas[iset][itrigg][iassoc]->GetBinContent(maxBin) + 10.*hCorrMeas[iset][itrigg][iassoc]->GetBinContent(maxBin);
                if (rangeMin<=0) rangeMin = 0.0008;
                cout << rangeMin << " " << rangeMax << endl;
                TPad *p = fpad[iset][itrigg][iassoc]->GetPad(1);
                p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
                hset(*hCorrReal[iset][itrigg][iassoc], "#Delta#phi", "1/N_{trigg}dN/d#Delta#phi", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
//                hset(*hCorrReal[iset][itrigg][iassoc], "#Delta#phi", "arb. norm.", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
                hCorrReal[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(rangeMin, rangeMax);
//                hCorrReal[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(TMath::Pi()/2., 3.*TMath::Pi()/2.);
                hCorrReal[iset][itrigg][iassoc]->Draw("P");
                hCorrMeas[iset][itrigg][iassoc]->Draw("SAME P");
                hCorrFinal[iset][itrigg][iassoc]->Draw("SAME P");
                leg[iset][itrigg][iassoc]->Draw("SAME");
                t->Draw("SAME");

                // Lower pad
                p = fpad[iset][itrigg][iassoc]->GetPad(2);
                p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
                hset( *hRatioMeas[iset][itrigg][iassoc], "#Delta#phi", "Ratio",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
                hRatioMeas[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.2, 1.8);
                //hRatioMeas[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(TMath::Pi()/2., 3.*TMath::Pi()/2.);
                hRatioMeas[iset][itrigg][iassoc]->Draw("P");
                hRatioFinal[iset][itrigg][iassoc]->Draw("P SAME");
                //hRatioSBMeas[iset][itrigg][iassoc]->Draw("P SAME");
            }
        }
    }

    //fpadPtComp = new Filipad(padID, 1.1, 0.35, 100, 100, 0.8, 1);
    //fpadPtComp->Draw();
    //padID++;

    //// Upper pad
    //TPad *p = fpadPtComp->GetPad(1);
    //p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
    //hset(*hCorrReal[0][0][0], "#Delta#phi", "1/N_{trigg}dN/d#phi", 1.1,1.5, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    //hCorrReal[0][0][0]->Draw("P");
    //hCorrReal[0][1][1]->Draw("SAME P");

    //// Lower pad
    //p = fpadPtComp->GetPad(2);
    //p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    //hset( *hRatioPt, "#Delta#phi", "Ratio",1.1,0.9, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
    //hRatioPt->Draw("P");
}
