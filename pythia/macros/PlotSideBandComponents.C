#include "include/Filipad.h"
#include "include/rootcommon.h"

void LoadData();
void ConfigHistos();
void DrawFiliPad();

const int nTriggBins = 2;
const int nAssocBins = 2;
//const double triggPt[nTriggBins+1] = {1.0, 2.0, 2.5, 3.0};
//const double assocPt[nAssocBins+1] = {0.5, 1.0, 1.5, 2.0, 2.5};
//const double triggPt[nTriggBins+1] = {4.0, 8.0, 8.0};
//const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};
const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

const double ntrigg[nTriggBins] = {39032, 3128}; // GEANT
//const double ntrigg[nTriggBins] = {249437, 10460}; // PYTHIA

const double mSize = 1.;

const int nset = 1;
//TString infiles[nset] = {
//    "analysis_STAR_pp.root"//,
//	//"analysis_STAR_pAu.root"
//};

TString infiles[nset] = {
    //"analysis_FoCal_pp_geant.root"
    "analysis_FoCal_pp_pythia.root"
    //"analysis_FoCal_pPb_pythia.root"
};

TString legHeader[nset] = {
    "p-p #sqrt{s} = 14 TeV"//,
    //"pPb #sqrt{s} = 5.02 TeV"//,
    //"p-Au"
};

TFile *fin[nset];
Filipad *fpad[nset][nTriggBins][nAssocBins];
Filipad *fpadPtComp;
TH1D *hCorrReal[nset][nTriggBins][nAssocBins];
TH1D *hCorrMeas[nset][nTriggBins][nAssocBins];
TH1D *hCorrMassMass[nset][nTriggBins][nAssocBins];
TH1D *hCorrMassSide[nset][nTriggBins][nAssocBins];
TH1D *hCorrSideMass[nset][nTriggBins][nAssocBins];
TH1D *hCorrSideSide[nset][nTriggBins][nAssocBins];
TH1D *hRatioMassMass[nset][nTriggBins][nAssocBins];
TH1D *hRatioMassSide[nset][nTriggBins][nAssocBins];
TH1D *hRatioSideMass[nset][nTriggBins][nAssocBins];
TH1D *hRatioSideSide[nset][nTriggBins][nAssocBins];
TLegend *leg[nset][nTriggBins][nAssocBins];


void PlotSideBandComponents()
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

                hCorrReal[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));  //hCorrReal[iset][itrigg][iassoc]->Rebin(6);
                hCorrMassMass[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));  //hCorrMassMass[iset][itrigg][iassoc]->Rebin(6);
                hCorrMassSide[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));  //hCorrMassSide[iset][itrigg][iassoc]->Rebin(6);
                hCorrSideMass[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));  //hCorrSideMass[iset][itrigg][iassoc]->Rebin(6);
                hCorrSideSide[iset][itrigg][iassoc]  = (TH1D*)fin[iset]->Get(Form("hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));  //hCorrSideSide[iset][itrigg][iassoc]->Rebin(6);

                cout << "bin width = " << hCorrMassMass[iset][itrigg][iassoc]->GetBinWidth(0) << endl;

                hCorrMeas[iset][itrigg][iassoc] = (TH1D*)hCorrMassMass[iset][itrigg][iassoc]->Clone(Form("hCorrMeas[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));

                hCorrMassMass[iset][itrigg][iassoc]->Scale(1./ntrigg[itrigg], "width");
                hCorrMassSide[iset][itrigg][iassoc]->Scale(1./ntrigg[itrigg], "width");
                hCorrSideMass[iset][itrigg][iassoc]->Scale(1./ntrigg[itrigg], "width");
                hCorrSideSide[iset][itrigg][iassoc]->Scale(1./ntrigg[itrigg], "width");

                // Ratios
                hRatioMassMass[iset][itrigg][iassoc] = (TH1D*)hCorrMassMass[iset][itrigg][iassoc]->Clone(Form("hRatioMassMass[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioMassMass[iset][itrigg][iassoc]->Divide(hCorrMassMass[iset][itrigg][iassoc]);
                hRatioMassSide[iset][itrigg][iassoc] = (TH1D*)hCorrMassSide[iset][itrigg][iassoc]->Clone(Form("hRatioMassSide[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioMassSide[iset][itrigg][iassoc]->Divide(hCorrMassMass[iset][itrigg][iassoc]);
                hRatioSideMass[iset][itrigg][iassoc] = (TH1D*)hCorrSideMass[iset][itrigg][iassoc]->Clone(Form("hRatioSideMass[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioSideMass[iset][itrigg][iassoc]->Divide(hCorrMassMass[iset][itrigg][iassoc]);
                hRatioSideSide[iset][itrigg][iassoc] = (TH1D*)hCorrSideSide[iset][itrigg][iassoc]->Clone(Form("hRatioSideSide[%4.1f,%4.1f][%4.1f,%4.1f]_px",tlow,tupp,alow,aupp));
                hRatioSideSide[iset][itrigg][iassoc]->Divide(hCorrMassMass[iset][itrigg][iassoc]);
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
                hCorrReal[iset][itrigg][iassoc]->SetMarkerStyle(kFullCircle);
                hCorrReal[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hCorrMeas[iset][itrigg][iassoc]->SetLineColor(kBlue);
                hCorrMeas[iset][itrigg][iassoc]->SetMarkerColor(kBlue);
                hCorrMeas[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hCorrMeas[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hCorrMassMass[iset][itrigg][iassoc]->SetLineColor(kBlack);
                hCorrMassMass[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrMassMass[iset][itrigg][iassoc]->SetMarkerStyle(kFullCircle);
                hCorrMassMass[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                //hCorrMassMass[iset][itrigg][iassoc]->Add(hCorrMassSide[iset][itrigg][iassoc], -1.87224);
                //hCorrMassMass[iset][itrigg][iassoc]->Add(hCorrSideMass[iset][itrigg][iassoc], -2.24933);
                //hCorrMassMass[iset][itrigg][iassoc]->Add(hCorrSideSide[iset][itrigg][iassoc], 4.21128);

                //hCorrMassMass[iset][itrigg][iassoc]->Add(hCorrMassSide[iset][itrigg][iassoc]);
                //hCorrMassMass[iset][itrigg][iassoc]->Add(hCorrSideMass[iset][itrigg][iassoc]);
                //hCorrMassMass[iset][itrigg][iassoc]->Add(hCorrSideSide[iset][itrigg][iassoc]);

                hCorrMassSide[iset][itrigg][iassoc]->SetLineColor(kRed);
                hCorrMassSide[iset][itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrMassSide[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hCorrMassSide[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hCorrSideMass[iset][itrigg][iassoc]->SetLineColor(kBlue);
                hCorrSideMass[iset][itrigg][iassoc]->SetMarkerColor(kBlue);
                hCorrSideMass[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hCorrSideMass[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hCorrSideSide[iset][itrigg][iassoc]->SetLineColor(kMagenta);
                hCorrSideSide[iset][itrigg][iassoc]->SetMarkerColor(kMagenta);
                hCorrSideSide[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hCorrSideSide[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hRatioMassMass[iset][itrigg][iassoc]->SetLineColor(kBlack);
                hRatioMassMass[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hRatioMassMass[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hRatioMassMass[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hRatioMassSide[iset][itrigg][iassoc]->SetLineColor(kRed);
                hRatioMassSide[iset][itrigg][iassoc]->SetMarkerColor(kRed);
                hRatioMassSide[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hRatioMassSide[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hRatioSideMass[iset][itrigg][iassoc]->SetLineColor(kBlue);
                hRatioSideMass[iset][itrigg][iassoc]->SetMarkerColor(kBlue);
                hRatioSideMass[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hRatioSideMass[iset][itrigg][iassoc]->SetMarkerSize(mSize);

                hRatioSideSide[iset][itrigg][iassoc]->SetLineColor(kMagenta);
                hRatioSideSide[iset][itrigg][iassoc]->SetMarkerColor(kMagenta);
                hRatioSideSide[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hRatioSideSide[iset][itrigg][iassoc]->SetMarkerSize(mSize);


                //leg[iset][itrigg][iassoc] = new TLegend(0.58, 0.55, 0.85, 0.78);
                leg[iset][itrigg][iassoc] = new TLegend(0.45, 0.6, 0.89, 0.82);
                //leg[iset][itrigg][iassoc] = new TLegend(0.25, 0.02, 0.75, 0.2);
                leg[iset][itrigg][iassoc]->SetNColumns(2);
                leg[iset][itrigg][iassoc]->SetFillStyle(0); leg[iset][itrigg][iassoc]->SetBorderSize(0); leg[iset][itrigg][iassoc]->SetTextSize(0.05);
                leg[iset][itrigg][iassoc]->SetHeader(Form("#splitline{%s}{[%0.1f,%0.1f][%0.1f,%0.1f]}", legHeader[iset].Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
                //leg[iset][itrigg][iassoc]->SetHeader(Form("%s, [%0.1f,%0.1f][%0.1f,%0.1f]", legHeader[iset].Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
                //leg[iset][itrigg][iassoc]->AddEntry(hCorrReal[0][0][0], "MC truth", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrMassMass[0][0][0], "mass-mass", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrMassSide[0][0][0], "mass-side", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrSideMass[0][0][0], "side-mass", "pe");
                leg[iset][itrigg][iassoc]->AddEntry(hCorrSideSide[0][0][0], "side-side", "pe");
            }
        }
    }
}

void DrawFiliPad()
{
    TText *t = new TText(.58,0.05,"PYTHIA8 simulation");
    //TText *t = new TText(0.58,0.78,"PYTHIA8 simulation");
    //TText *t = new TText(.58,0.78,"PYTHIA6+GEANT3");
    t->SetNDC();
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
                int minBin = hCorrMassMass[iset][itrigg][iassoc]->GetMinimumBin();
                int maxBin = hCorrMassMass[iset][itrigg][iassoc]->GetMaximumBin();
                double rangeMin = hCorrMassMass[iset][itrigg][iassoc]->GetBinContent(minBin) - 0.4*hCorrMassMass[iset][itrigg][iassoc]->GetBinContent(minBin);
                double rangeMax = hCorrMassMass[iset][itrigg][iassoc]->GetBinContent(maxBin) + 2.*hCorrMassMass[iset][itrigg][iassoc]->GetBinContent(maxBin);
                if (rangeMin<=0) rangeMin = 0.0008;
                rangeMin = 5*10e-6;
                rangeMax = 20;
                cout << rangeMin << " " << rangeMax << endl;
                TPad *p = fpad[iset][itrigg][iassoc]->GetPad(1);
                p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
//                hset(*hCorrReal[iset][itrigg][iassoc], "#Delta#phi", "1/N_{trigg}dN/d#phi", 1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
                hset(*hCorrMassMass[iset][itrigg][iassoc], "#Delta#phi", "1/N_{trigg}dN/d#Delta#phi", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
                //hCorrMassMass[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.1, 0.1);
                hCorrMassMass[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(rangeMin, rangeMax);
                hCorrMassMass[iset][itrigg][iassoc]->Draw("P");
                //hCorrReal[iset][itrigg][iassoc]->Draw("SAME P");
                //hCorrMeas[iset][itrigg][iassoc]->Draw("SAME P");
                hCorrMassSide[iset][itrigg][iassoc]->Draw("SAME P");
                hCorrSideMass[iset][itrigg][iassoc]->Draw("SAME P");
                hCorrSideSide[iset][itrigg][iassoc]->Draw("SAME P");
                leg[iset][itrigg][iassoc]->Draw("SAME");
                t->Draw("SAME");

                // Lower pad
                p = fpad[iset][itrigg][iassoc]->GetPad(2);
                p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
                hset( *hRatioMassSide[iset][itrigg][iassoc], "#Delta#phi", "Ratio",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
                //hRatioMassMass[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.2, 1.8);
                //hRatioMassMass[iset][itrigg][iassoc]->Draw("P");
                hRatioMassSide[iset][itrigg][iassoc]->Draw("P");
                hRatioSideMass[iset][itrigg][iassoc]->Draw("SAME P");
                hRatioSideSide[iset][itrigg][iassoc]->Draw("SAME P");
            }
        }
    }
}
