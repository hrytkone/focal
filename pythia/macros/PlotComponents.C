#include "include/Filipad.h"
#include "include/rootcommon.h"

const bool plotCloseup = 0;
const bool useLeading = 0;

const int nTriggBins = 2;
const int nAssocBins = 2;
const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

TString legHeader = "p-p #sqrt{s} = 14 TeV";
int nEvent = 1;

TFile *fIn;
TH1D *hCounter;
Filipad *fpad[nTriggBins][nAssocBins];
Filipad *fpadcomp[nTriggBins][nAssocBins];

TH2D* hCorrTrue[nTriggBins][nAssocBins];
TH1D* hCorrTrueProj[nTriggBins][nAssocBins];

TH2D* hCorrMassMass[nTriggBins][nAssocBins];
TH2D* hCorrMassSide[nTriggBins][nAssocBins];
TH2D* hCorrSideMass[nTriggBins][nAssocBins];
TH2D* hCorrSideSide[nTriggBins][nAssocBins];

TH1D* hCorrMassMassProj[nTriggBins][nAssocBins];
TH1D* hCorrMassSideProj[nTriggBins][nAssocBins];
TH1D* hCorrSideMassProj[nTriggBins][nAssocBins];
TH1D* hCorrSideSideProj[nTriggBins][nAssocBins];

TH2D* hCorrSS[nTriggBins][nAssocBins];
TH2D* hCorrSB[nTriggBins][nAssocBins];
TH2D* hCorrBS[nTriggBins][nAssocBins];
TH2D* hCorrBB[nTriggBins][nAssocBins];

TH1D* hCorrSSProj[nTriggBins][nAssocBins];
TH1D* hCorrSBProj[nTriggBins][nAssocBins];
TH1D* hCorrBSProj[nTriggBins][nAssocBins];
TH1D* hCorrBBProj[nTriggBins][nAssocBins];

TH1D* hCorrSSRec[nTriggBins][nAssocBins];
TH1D* hCorrSBRec[nTriggBins][nAssocBins];
TH1D* hCorrBSRec[nTriggBins][nAssocBins];
TH1D* hCorrBBRec[nTriggBins][nAssocBins];

TH1D *hTrueVsMeasured[nTriggBins][nAssocBins];
TH1D *hSSVsTrue[nTriggBins][nAssocBins];
TH1D *hSSVsMeasured[nTriggBins][nAssocBins];
TH1D *hSBVsMeasured[nTriggBins][nAssocBins];
TH1D *hBSVsMeasured[nTriggBins][nAssocBins];
TH1D *hBBVsMeasured[nTriggBins][nAssocBins];

TLegend *leg1[nTriggBins][nAssocBins];
TLegend *leg2[nTriggBins][nAssocBins];
TLegend *leg3[nTriggBins][nAssocBins];

void LoadData(TString input);
void ConfigHistos();
void DrawFiliPad();
void CreateLegends();

void PlotComponents(TString sInputName = "input.root")
{
    gStyle->SetOptStat(0);
    LoadData(sInputName);
    ConfigHistos();
    CreateLegends();
    DrawFiliPad();
}


// -----------------
// |   Functions   |
// -----------------

void LoadData(TString input)
{
    fIn = TFile::Open(input);
    hCounter = (TH1D*)fIn->Get("hCounter");
    nEvent = hCounter->GetBinContent(1);
    std::cout << "Input file contains " << nEvent << " events, proceed to analyse" << std::endl;

    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (!useLeading && tlow < aupp) continue;
            hCorrTrue[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrTrue[itrigg][iassoc]->Rebin2D(4);
            hCorrTrue[itrigg][iassoc]->Scale(1./nEvent);
            hCorrTrue[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrMassMass[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassMass[itrigg][iassoc]->Rebin2D(4);
            hCorrMassMass[itrigg][iassoc]->Scale(1./nEvent);
            hCorrMassMass[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrMassSide[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassSide[itrigg][iassoc]->Rebin2D(4);
            hCorrMassSide[itrigg][iassoc]->Scale(1./nEvent);
            hCorrMassSide[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSideMass[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideMass[itrigg][iassoc]->Rebin2D(4);
            hCorrSideMass[itrigg][iassoc]->Scale(1./nEvent);
            hCorrSideMass[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSideSide[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideSide[itrigg][iassoc]->Rebin2D(4);
            hCorrSideSide[itrigg][iassoc]->Scale(1./nEvent);
            hCorrSideSide[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSS[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrSignalSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSS[itrigg][iassoc]->Rebin2D(4);
            hCorrSS[itrigg][iassoc]->Scale(1./nEvent);
            hCorrSS[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSB[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrSignalBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSB[itrigg][iassoc]->Rebin2D(4);
            hCorrSB[itrigg][iassoc]->Scale(1./nEvent);
            hCorrSB[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrBS[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrBgSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrBS[itrigg][iassoc]->Rebin2D(4);
            hCorrBS[itrigg][iassoc]->Scale(1./nEvent);
            hCorrBS[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrBB[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrBgBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrBB[itrigg][iassoc]->Rebin2D(4);
            hCorrBB[itrigg][iassoc]->Scale(1./nEvent);
            hCorrBB[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);
        }
    }
}

void ConfigHistos()
{
    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (!useLeading && tlow < aupp) continue;
                hCorrTrueProj[itrigg][iassoc] = hCorrTrue[itrigg][iassoc]->ProjectionX();
                hCorrTrueProj[itrigg][iassoc]->SetTitle("; #Delta#phi; 1/N dN/d#Delta#phi");
                hCorrTrueProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                //hCorrTrueProj[itrigg][iassoc]->GetYaxis()->SetRangeUser(0.00001, 0.6);
                hCorrTrueProj[itrigg][iassoc]->SetLineColor(kBlack);
                hCorrTrueProj[itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrTrueProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrTrueProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrTrueProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrMassMassProj[itrigg][iassoc] = hCorrMassMass[itrigg][iassoc]->ProjectionX();
                hCorrMassMassProj[itrigg][iassoc]->SetTitle("; #Delta#phi; 1/N dN/d#Delta#phi");
                hCorrMassMassProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                //hCorrMassMassProj[itrigg][iassoc]->GetYaxis()->SetRangeUser(0.00001, 0.6);
                hCorrMassMassProj[itrigg][iassoc]->SetLineColor(kRed);
                hCorrMassMassProj[itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrMassMassProj[itrigg][iassoc]->SetMarkerStyle(20);
                hCorrMassMassProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrMassMassProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrMassSideProj[itrigg][iassoc] = hCorrMassSide[itrigg][iassoc]->ProjectionX();
                hCorrMassSideProj[itrigg][iassoc]->SetTitle("; #Delta#phi; 1/N dN/d#Delta#phi");
                hCorrMassSideProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hCorrMassSideProj[itrigg][iassoc]->SetLineColor(kRed);
                hCorrMassSideProj[itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrMassSideProj[itrigg][iassoc]->SetMarkerStyle(20);
                hCorrMassSideProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrMassSideProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrSideMassProj[itrigg][iassoc] = hCorrSideMass[itrigg][iassoc]->ProjectionX();
                hCorrSideMassProj[itrigg][iassoc]->SetTitle("; #Delta#phi; 1/N dN/d#Delta#phi");
                hCorrSideMassProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hCorrSideMassProj[itrigg][iassoc]->SetLineColor(kRed);
                hCorrSideMassProj[itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrSideMassProj[itrigg][iassoc]->SetMarkerStyle(20);
                hCorrSideMassProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrSideMassProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrSideSideProj[itrigg][iassoc] = hCorrSideSide[itrigg][iassoc]->ProjectionX();
                hCorrSideSideProj[itrigg][iassoc]->SetTitle("; #Delta#phi; 1/N dN/d#Delta#phi");
                hCorrSideSideProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hCorrSideSideProj[itrigg][iassoc]->SetLineColor(kRed);
                hCorrSideSideProj[itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrSideSideProj[itrigg][iassoc]->SetMarkerStyle(20);
                hCorrSideSideProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrSideSideProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrSSProj[itrigg][iassoc] = hCorrSS[itrigg][iassoc]->ProjectionX();
                hCorrSSProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrSSProj[itrigg][iassoc]->SetMarkerColor(kMagenta-7);
                hCorrSSProj[itrigg][iassoc]->SetLineColor(kMagenta-7);
                hCorrSSProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrSSProj[itrigg][iassoc]->SetLineWidth(2);
                hCorrSSProj[itrigg][iassoc]->SetTitle("; #Delta#phi; 1/N dN/d#Delta#phi");

                hCorrSBProj[itrigg][iassoc] = hCorrSB[itrigg][iassoc]->ProjectionX();
                hCorrSBProj[itrigg][iassoc]->SetMarkerStyle(7);
                hCorrSBProj[itrigg][iassoc]->SetMarkerColor(kOrange+7);
                hCorrSBProj[itrigg][iassoc]->SetLineColor(kOrange+7);
                hCorrSBProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrSBProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrBSProj[itrigg][iassoc] = hCorrBS[itrigg][iassoc]->ProjectionX();
                hCorrBSProj[itrigg][iassoc]->SetMarkerStyle(7);
                hCorrBSProj[itrigg][iassoc]->SetMarkerColor(kCyan+1);
                hCorrBSProj[itrigg][iassoc]->SetLineColor(kCyan+1);
                hCorrBSProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrBSProj[itrigg][iassoc]->SetTitle("; #Delta#phi; 1/N dN/d#Delta#phi");
                hCorrBSProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrBBProj[itrigg][iassoc] = hCorrBB[itrigg][iassoc]->ProjectionX();
                hCorrBBProj[itrigg][iassoc]->SetMarkerStyle(7);
                hCorrBBProj[itrigg][iassoc]->SetMarkerColor(kBlue-7);
                hCorrBBProj[itrigg][iassoc]->SetLineColor(kBlue-7);
                hCorrBBProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrBBProj[itrigg][iassoc]->SetLineWidth(2);

                // RATIOS
                hTrueVsMeasured[itrigg][iassoc] = (TH1D*)hCorrTrueProj[itrigg][iassoc]->Clone(Form("hTrueVsMeasured[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hTrueVsMeasured[itrigg][iassoc]->Divide(hCorrMassMassProj[itrigg][iassoc]);
                hSSVsMeasured[itrigg][iassoc] = (TH1D*)hCorrSSProj[itrigg][iassoc]->Clone(Form("hSSVsMeasured[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hSSVsMeasured[itrigg][iassoc]->Divide(hCorrMassMassProj[itrigg][iassoc]);
                hSSVsTrue[itrigg][iassoc] = (TH1D*)hCorrSSProj[itrigg][iassoc]->Clone(Form("hSSVsTrue[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hSSVsTrue[itrigg][iassoc]->Divide(hCorrTrueProj[itrigg][iassoc]);
                hSSVsTrue[itrigg][iassoc]->SetLineColor(kMagenta-3);

                hSBVsMeasured[itrigg][iassoc] = (TH1D*)hCorrSBProj[itrigg][iassoc]->Clone(Form("hSBVsMeasured[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hSBVsMeasured[itrigg][iassoc]->Divide(hCorrMassMassProj[itrigg][iassoc]);
                hBSVsMeasured[itrigg][iassoc] = (TH1D*)hCorrBSProj[itrigg][iassoc]->Clone(Form("hSBVsMeasured[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hBSVsMeasured[itrigg][iassoc]->Divide(hCorrMassMassProj[itrigg][iassoc]);
                hBBVsMeasured[itrigg][iassoc] = (TH1D*)hCorrBBProj[itrigg][iassoc]->Clone(Form("hSBVsMeasured[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hBBVsMeasured[itrigg][iassoc]->Divide(hCorrMassMassProj[itrigg][iassoc]);
            }
        }
}

void CreateLegends()
{
    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (!useLeading && tlow < aupp) continue;

            if (plotCloseup) {
                leg1[itrigg][iassoc] = new TLegend(0.25, 0.65, 0.85, 0.78);
                leg1[itrigg][iassoc]->SetNColumns(3);
                leg1[itrigg][iassoc]->SetFillStyle(0); leg1[itrigg][iassoc]->SetBorderSize(0); leg1[itrigg][iassoc]->SetTextSize(0.05);
                leg1[itrigg][iassoc]->SetHeader(Form("%s, [%0.1f,%0.1f][%0.1f,%0.1f]", legHeader.Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));

                leg2[itrigg][iassoc] = new TLegend(0.25, 0.3, 0.9, 0.55);
                leg2[itrigg][iassoc]->SetNColumns(2);
                leg2[itrigg][iassoc]->SetFillStyle(0); leg2[itrigg][iassoc]->SetBorderSize(0); leg2[itrigg][iassoc]->SetTextSize(0.06);

                leg3[itrigg][iassoc] = new TLegend(0.25, 0.6, 0.85, 0.78);
                leg3[itrigg][iassoc]->SetFillStyle(0); leg3[itrigg][iassoc]->SetBorderSize(0); leg3[itrigg][iassoc]->SetTextSize(0.05);
                leg3[itrigg][iassoc]->SetNColumns(3);
                leg3[itrigg][iassoc]->SetHeader(Form("%s, [%0.1f,%0.1f][%0.1f,%0.1f]", legHeader.Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
            } else {
                leg1[itrigg][iassoc] = new TLegend(0.5, 0.4, 0.78, 0.75);
                leg1[itrigg][iassoc]->SetFillStyle(0); leg1[itrigg][iassoc]->SetBorderSize(0); leg1[itrigg][iassoc]->SetTextSize(0.05);
                leg1[itrigg][iassoc]->SetHeader(Form("#splitline{%s}{[%0.1f,%0.1f][%0.1f,%0.1f]}", legHeader.Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));

                leg2[itrigg][iassoc] = new TLegend(0.6, 0.3, 0.9, 0.6);
                leg2[itrigg][iassoc]->SetFillStyle(0); leg2[itrigg][iassoc]->SetBorderSize(0); leg2[itrigg][iassoc]->SetTextSize(0.06);

                leg3[itrigg][iassoc] = new TLegend(0.48, 0.5, 0.92, 0.75);
                leg3[itrigg][iassoc]->SetFillStyle(0); leg3[itrigg][iassoc]->SetBorderSize(0); leg3[itrigg][iassoc]->SetTextSize(0.05);
                leg3[itrigg][iassoc]->SetNColumns(3);
                leg3[itrigg][iassoc]->SetHeader(Form("#splitline{%s}{[%0.1f,%0.1f][%0.1f,%0.1f]}", legHeader.Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
            }
            leg1[itrigg][iassoc]->AddEntry(hCorrTrueProj[itrigg][iassoc], "MC truth", "le");
            leg1[itrigg][iassoc]->AddEntry(hCorrMassMassProj[itrigg][iassoc], "f_{mass,mass}", "le");
            leg1[itrigg][iassoc]->AddEntry(hCorrSSProj[itrigg][iassoc], "f_{true,true}", "le");

            leg2[itrigg][iassoc]->AddEntry(hTrueVsMeasured[itrigg][iassoc], "MC truth / f(mass,mass)", "le");
            leg2[itrigg][iassoc]->AddEntry(hSSVsMeasured[itrigg][iassoc], "f(true,true) / f(mass,mass)", "le");
            leg2[itrigg][iassoc]->AddEntry(hSSVsTrue[itrigg][iassoc], "f(true,true) / MC truth", "le");

            leg3[itrigg][iassoc]->AddEntry(hCorrMassMassProj[itrigg][iassoc], "f_{mass,mass}", "le");
            leg3[itrigg][iassoc]->AddEntry(hCorrSSProj[itrigg][iassoc], "f_{true,true}", "le");
            leg3[itrigg][iassoc]->AddEntry(hCorrSBProj[itrigg][iassoc], "f_{true,fake}", "le");
            leg3[itrigg][iassoc]->AddEntry(hCorrBSProj[itrigg][iassoc], "f_{fake,true}", "le");
            leg3[itrigg][iassoc]->AddEntry(hCorrBBProj[itrigg][iassoc], "f_{fake,fake}", "le");
        }
    }
}

void DrawFiliPad()
{
    TText *t;
    if (plotCloseup)
        t = new TText(.25,0.79,"PYTHIA8 simulation");
    else
        t = new TText(.48,0.79,"PYTHIA8 simulation");
    t->SetNDC();
    int padID = 0;
    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (!useLeading && tlow < aupp) continue;

//                fpad[itrigg][iassoc] = new Filipad(padID, 1.1, 0.5, 100, 100, 0.7, 1);
            fpad[itrigg][iassoc] = new Filipad(padID, 1.1, 0.35, 100, 100, 0.8, 1);
            fpad[itrigg][iassoc]->Draw();
            padID++;

            // Upper pad
            int minBin = hCorrSSProj[itrigg][iassoc]->GetMinimumBin();
            int maxBin = hCorrMassMassProj[itrigg][iassoc]->GetMaximumBin();
            double rangeMin = hCorrSSProj[itrigg][iassoc]->GetBinContent(minBin) - 0.2*hCorrSSProj[itrigg][iassoc]->GetBinContent(minBin);
            double rangeMax = hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin) + hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin);
            if (rangeMin<=0) rangeMin = 0.0008;
            if (itrigg==1 && iassoc==1) {
                rangeMin = 0.5;
                rangeMax *= 2.;
            }
            if (plotCloseup) {
                rangeMax = hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin) + 10.*hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin);
            }
            rangeMin = 5e-7;
            rangeMax = 2e-4;
            TPad *p = fpad[itrigg][iassoc]->GetPad(1);
            p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
            hset(*hCorrTrueProj[itrigg][iassoc], "#Delta#phi", "1/N_{ev} dN/d#Delta#phi", 1.1,1.4, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
            hCorrTrueProj[itrigg][iassoc]->GetYaxis()->SetRangeUser(rangeMin, rangeMax);
            hCorrTrueProj[itrigg][iassoc]->Draw("HIST E");
            hCorrMassMassProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrSSProj[itrigg][iassoc]->Draw("SAME HIST E");
            leg1[itrigg][iassoc]->Draw("SAME");
            t->Draw("SAME");

            // Lower pad
            p = fpad[itrigg][iassoc]->GetPad(2);
            p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
            hset( *hTrueVsMeasured[itrigg][iassoc], "#Delta#phi", "Ratio",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
            if (plotCloseup)
                hTrueVsMeasured[itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.8, 1.1);
            else
                hTrueVsMeasured[itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.8, 1.8);
            hTrueVsMeasured[itrigg][iassoc]->Draw("HIST E");
            hSSVsMeasured[itrigg][iassoc]->Draw("HIST E SAME");
            hSSVsTrue[itrigg][iassoc]->Draw("HIST E SAME");
            leg2[itrigg][iassoc]->Draw("SAME");
        }
    }

    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (!useLeading && tlow < aupp) continue;

//                fpad[itrigg][iassoc] = new Filipad(padID, 1.1, 0.5, 100, 100, 0.7, 1);
            fpadcomp[itrigg][iassoc] = new Filipad(padID, 1.1, 0.35, 100, 100, 0.8, 1);
            fpadcomp[itrigg][iassoc]->Draw();
            padID++;

            // Upper pad
            int minBin = hCorrSSProj[itrigg][iassoc]->GetMinimumBin();
            int maxBin = hCorrMassMassProj[itrigg][iassoc]->GetMaximumBin();
            double rangeMin = 0.5;
            double rangeMax = hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin) + 10.*hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin);
            if (rangeMin<=0) rangeMin = 0.0000008;
            if (itrigg==0 && iassoc==0) {
                rangeMin = 1;
                rangeMax *= 3.;
            }
            if (plotCloseup) {
                rangeMax = hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin) + 100.*hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin);
            }
            if (iassoc==0) {
                rangeMin = 5e-9;
            } else {
                rangeMin = 5e-10;
            }
            rangeMax = 5e-4;
            TPad *p = fpadcomp[itrigg][iassoc]->GetPad(1);
            p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
            hset(*hCorrMassMassProj[itrigg][iassoc], "#Delta#phi", "1/N_{ev} dN/d#Delta#phi", 1.2,1.4, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
            hCorrMassMassProj[itrigg][iassoc]->GetYaxis()->SetRangeUser(rangeMin, rangeMax);
            hCorrMassMassProj[itrigg][iassoc]->Draw("HIST E");
            hCorrSSProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrSBProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrBSProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrBBProj[itrigg][iassoc]->Draw("SAME HIST E");
            leg3[itrigg][iassoc]->Draw("SAME");
            t->Draw("SAME");

            // Lower pad
            p = fpadcomp[itrigg][iassoc]->GetPad(2);
            p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
            hset( *hSSVsMeasured[itrigg][iassoc], "#Delta#phi", "f_{..,..}/f_{mass,mass}",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
            hSSVsMeasured[itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.2, .95);
            hSSVsMeasured[itrigg][iassoc]->Draw("HIST E");
            hSBVsMeasured[itrigg][iassoc]->Draw("SAME HIST E");
            hBSVsMeasured[itrigg][iassoc]->Draw("SAME HIST E");
            hBBVsMeasured[itrigg][iassoc]->Draw("SAME HIST E");
        }
    }

    /**for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (!useLeading && tlow < aupp) continue;

            TCanvas *c1 = new TCanvas(Form("c_%d-%d", itrigg, iassoc), "", 800, 800);
            c1->Divide(2,2);

            c1->cd(1);
            hCorrSSProj[itrigg][iassoc]->Draw("HIST E");

            c1->cd(2);
            hCorrMassSideProj[itrigg][iassoc]->Draw("HIST E");
            hCorrMassSideProj[itrigg][iassoc]->Scale(1./hCorrMassSideProj[itrigg][iassoc]->GetEntries());
            hCorrSBProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrSBProj[itrigg][iassoc]->Scale(1./hCorrSBProj[itrigg][iassoc]->GetEntries());
            hCorrBBProj[itrigg][iassoc]->Draw("SAME HIST E");

            c1->cd(3);
            hCorrSideMassProj[itrigg][iassoc]->Draw("HIST E");
            hCorrSideMassProj[itrigg][iassoc]->Scale(1./hCorrSideSideProj[itrigg][iassoc]->GetEntries());
            hCorrBSProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrBSProj[itrigg][iassoc]->Scale(1./hCorrBSProj[itrigg][iassoc]->GetEntries());
            hCorrBBProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrBBProj[itrigg][iassoc]->Scale(1./hCorrBBProj[itrigg][iassoc]->GetEntries());

            c1->cd(4);
            hCorrSideSideProj[itrigg][iassoc]->Draw("HIST E");
            hCorrSideSideProj[itrigg][iassoc]->Scale(1./hCorrSideSideProj[itrigg][iassoc]->GetEntries());
            hCorrBBProj[itrigg][iassoc]->Draw("SAME HIST E");
        }
    }**/
}
