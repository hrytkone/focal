#include "include/Filipad.h"
#include "include/rootcommon.h"

void LoadData();
void ConfigHistos();
void DrawFiliPad();

const int nTriggBins = 2;
const int nAssocBins = 2;
//double triggPt[nTriggBins+1] = {2.0, 2.5, 3.0};
//double assocPt[nAssocBins+1] = {1.0, 1.5, 2.0, 2.5};
const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

const int nset = 1;
//TString infiles[nset] = {
//    "analysis_STAR_pp.root"//,
//	//"analysis_STAR_pAu.root"
//};

TString infiles[nset] = {
    "analysis_FoCal_pp.root"
};

TString legHeader[nset] = {
    "p-p"//,
    //"p-Au"
};

TFile *fin[nset];
Filipad *fpad[nset][nTriggBins][nAssocBins];
TH1D *hCorrReal[nset][nTriggBins][nAssocBins];
TH1D *hCorrFinal[nset][nTriggBins][nAssocBins];
TH1D *hCorrMeas[nset][nTriggBins][nAssocBins];
TH1D *hRatio[nset][nTriggBins][nAssocBins];
TLegend *leg[nset][nTriggBins][nAssocBins];

void ComparisonPlot()
{
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

                hCorrReal[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp)); hCorrReal[iset][itrigg][iassoc]->Rebin(2);
                hCorrFinal[iset][itrigg][iassoc] = (TH1D*)fin[iset]->Get(Form("hCorrFinal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrFinal[iset][itrigg][iassoc]->Rebin(2);
				hCorrMeas[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrMeas[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp)); hCorrMeas[iset][itrigg][iassoc]->Rebin(2);

                // Calculate ratios
                hRatio[iset][itrigg][iassoc] = (TH1D*)hCorrReal[iset][itrigg][iassoc]->Clone(Form("hRatio[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatio[iset][itrigg][iassoc]->Divide(hCorrFinal[iset][itrigg][iassoc]);
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

                hCorrReal[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerSize(0.5);
                hCorrFinal[iset][itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrFinal[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hCorrFinal[iset][itrigg][iassoc]->SetMarkerSize(0.5);
                hCorrMeas[iset][itrigg][iassoc]->SetMarkerColor(kBlue);
                hCorrMeas[iset][itrigg][iassoc]->SetMarkerStyle(20);
                hCorrMeas[iset][itrigg][iassoc]->SetMarkerSize(0.5);
                hRatio[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(0., 5.);

                leg[iset][itrigg][iassoc] = new TLegend(0.68, 0.75, 0.88, 0.85);
                leg[iset][itrigg][iassoc]->SetFillStyle(0); leg[iset][itrigg][iassoc]->SetBorderSize(0); leg[iset][itrigg][iassoc]->SetTextSize(0.05);
                leg[iset][itrigg][iassoc]->SetHeader(Form("%s, [%0.1f,%0.1f][%0.1f,%0.1f]", legHeader[iset].Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
                leg[iset][itrigg][iassoc]->AddEntry(hCorrReal[0][0][0], "Real", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrMeas[0][0][0], "Meas", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrFinal[0][0][0], "Recovered", "pe");
            }
        }
    }
}

void DrawFiliPad()
{
    int padID = 0;
    for (int iset = 0; iset < nset; iset++) {
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;

                fpad[iset][itrigg][iassoc] = new Filipad(padID, 1.1, 0.5, 100, 100, 0.7, 1);
                fpad[iset][itrigg][iassoc]->Draw();
                padID++;

                // Upper pad
                TPad *p = fpad[iset][itrigg][iassoc]->GetPad(1);
                p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
                hset(*hCorrReal[iset][itrigg][iassoc], "#Delta#phi", "1/N_{trigg}dN/d#phi", 1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
                hCorrMeas[iset][itrigg][iassoc]->Draw("P");
                hCorrReal[iset][itrigg][iassoc]->Draw("SAME P");
                hCorrFinal[iset][itrigg][iassoc]->Draw("SAME P");
                leg[iset][itrigg][iassoc]->Draw("SAME");

                // Lower pad
                p = fpad[iset][itrigg][iassoc]->GetPad(2);
                p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
                hset( *hRatio[iset][itrigg][iassoc], "#Delta#phi", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
                hRatio[iset][itrigg][iassoc]->Draw();
            }
        }
    }
}
