#include "include/Filipad.h"
#include "include/rootcommon.h"

const bool plotCloseup = 0;
const bool useLeading = 0;

const int nTriggBins = 2;
const int nAssocBins = 2;
const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

TString legHeader = "p-p #sqrt{s} = 14 TeV";

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

void CheckBgComponents(TString sInputName = "input.root")
{
    gStyle->SetOptStat(0);
    LoadData(sInputName);
    ConfigHistos();
    DrawFiliPad();
}


// -----------------
// |   Functions   |
// -----------------

void LoadData(TString input)
{
    fIn = TFile::Open(input);
    hCounter = (TH1D*)fIn->Get("hCounter");
    int nEvent = hCounter->GetBinContent(1);
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
            hCorrTrue[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrMassMass[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassMass[itrigg][iassoc]->Rebin2D(4);
            hCorrMassMass[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrMassSide[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassSide[itrigg][iassoc]->Rebin2D(4);
            hCorrMassSide[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSideMass[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideMass[itrigg][iassoc]->Rebin2D(4);
            hCorrSideMass[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSideSide[itrigg][iassoc] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideSide[itrigg][iassoc]->Rebin2D(4);
            hCorrSideSide[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSS[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrSignalSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSS[itrigg][iassoc]->Rebin2D(4);
            hCorrSS[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrSB[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrSignalBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSB[itrigg][iassoc]->Rebin2D(4);
            hCorrSB[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrBS[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrBgSignal[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrBS[itrigg][iassoc]->Rebin2D(4);
            hCorrBS[itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);

            hCorrBB[itrigg][iassoc] = (TH2D*)fIn->Get(Form("TrueComponents/hCorrBgBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrBB[itrigg][iassoc]->Rebin2D(4);
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
                //hCorrTrueProj[itrigg][iassoc]->Scale(1./nMeasTrigg);
                hCorrTrueProj[itrigg][iassoc]->SetTitle("; #Delta#phi; Counts");
                hCorrTrueProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                //hCorrTrueProj[itrigg][iassoc]->GetYaxis()->SetRangeUser(0.00001, 0.6);
                hCorrTrueProj[itrigg][iassoc]->SetLineColor(kBlack);
                hCorrTrueProj[itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrTrueProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrTrueProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrTrueProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrMassMassProj[itrigg][iassoc] = hCorrMassMass[itrigg][iassoc]->ProjectionX();
                //hCorrMassMassProj[itrigg][iassoc]->Scale(1./nMeasTrigg);
                hCorrMassMassProj[itrigg][iassoc]->SetTitle("; #Delta#phi; Counts");
                hCorrMassMassProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                //hCorrMassMassProj[itrigg][iassoc]->GetYaxis()->SetRangeUser(0.00001, 0.6);
                hCorrMassMassProj[itrigg][iassoc]->SetLineColor(kRed);
                hCorrMassMassProj[itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrMassMassProj[itrigg][iassoc]->SetMarkerStyle(20);
                hCorrMassMassProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrMassMassProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrMassSideProj[itrigg][iassoc] = hCorrMassSide[itrigg][iassoc]->ProjectionX();
                hCorrMassSideProj[itrigg][iassoc]->SetTitle("; #Delta#phi; Counts");
                hCorrMassSideProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hCorrMassSideProj[itrigg][iassoc]->SetLineColor(kRed);
                hCorrMassSideProj[itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrMassSideProj[itrigg][iassoc]->SetMarkerStyle(20);
                hCorrMassSideProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrMassSideProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrSideMassProj[itrigg][iassoc] = hCorrSideMass[itrigg][iassoc]->ProjectionX();
                hCorrSideMassProj[itrigg][iassoc]->SetTitle("; #Delta#phi; Counts");
                hCorrSideMassProj[itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hCorrSideMassProj[itrigg][iassoc]->SetLineColor(kRed);
                hCorrSideMassProj[itrigg][iassoc]->SetMarkerColor(kRed);
                hCorrSideMassProj[itrigg][iassoc]->SetMarkerStyle(20);
                hCorrSideMassProj[itrigg][iassoc]->SetMarkerSize(0.);
                hCorrSideMassProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrSideSideProj[itrigg][iassoc] = hCorrSideSide[itrigg][iassoc]->ProjectionX();
                hCorrSideSideProj[itrigg][iassoc]->SetTitle("; #Delta#phi; Counts");
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
                //hCorrSSProj[itrigg][iassoc]->Scale(1./nRealRecTrigg);
                hCorrSSProj[itrigg][iassoc]->SetTitle("; #Delta#phi; Counts");

                hCorrSBProj[itrigg][iassoc] = hCorrSB[itrigg][iassoc]->ProjectionX();
                hCorrSBProj[itrigg][iassoc]->SetMarkerStyle(7);
                //hCorrSBProj[itrigg][iassoc]->Scale(1./nRealRecTrigg);
                hCorrSBProj[itrigg][iassoc]->SetMarkerColor(kOrange+7);
                hCorrSBProj[itrigg][iassoc]->SetLineColor(kOrange+7);
                hCorrSBProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrSBProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrBSProj[itrigg][iassoc] = hCorrBS[itrigg][iassoc]->ProjectionX();
                hCorrBSProj[itrigg][iassoc]->SetMarkerStyle(7);
                //hCorrBSProj[itrigg][iassoc]->Scale(1./nFakeRecTrigg);
                hCorrBSProj[itrigg][iassoc]->SetMarkerColor(kCyan+1);
                hCorrBSProj[itrigg][iassoc]->SetLineColor(kCyan+1);
                hCorrBSProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrBSProj[itrigg][iassoc]->SetTitle("; #Delta#phi; Counts");
                hCorrBSProj[itrigg][iassoc]->SetLineWidth(2);

                hCorrBBProj[itrigg][iassoc] = hCorrBB[itrigg][iassoc]->ProjectionX();
                hCorrBBProj[itrigg][iassoc]->SetMarkerStyle(7);
                //hCorrBBProj[itrigg][iassoc]->Scale(1./nFakeRecTrigg);
                hCorrBBProj[itrigg][iassoc]->SetMarkerColor(kBlue-7);
                hCorrBBProj[itrigg][iassoc]->SetLineColor(kBlue-7);
                hCorrBBProj[itrigg][iassoc]->SetMarkerStyle(kDot);
                hCorrBBProj[itrigg][iassoc]->SetLineWidth(2);
            }
        }
}

void DrawFiliPad()
{
    TText *t;
    if (plotCloseup)
        t = new TText(.25,0.79,"PYTHIA8 simulation");
    else
        t = new TText(.5,0.79,"PYTHIA8 simulation");
    t->SetNDC();

    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (!useLeading && tlow < aupp) continue;

            TCanvas *c1 = new TCanvas(Form("c_%d-%d", itrigg, iassoc), "", 800, 800);
            c1->Divide(2,2,0.0001,0.0001);

            c1->cd(1);
            hCorrSSProj[itrigg][iassoc]->Draw("HIST E");

            c1->cd(2);
            hCorrMassSideProj[itrigg][iassoc]->Draw("HIST E");
            hCorrSBProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrBBProj[itrigg][iassoc]->Draw("SAME HIST E");

            c1->cd(3);
            hCorrSideMassProj[itrigg][iassoc]->Draw("HIST E");
            hCorrBSProj[itrigg][iassoc]->Draw("SAME HIST E");
            hCorrBBProj[itrigg][iassoc]->Draw("SAME HIST E");

            c1->cd(4);
            //gPad->SetLogy();
            //hCorrSideSideProj[itrigg][iassoc]->Scale(1./hCorrSideSideProj[itrigg][iassoc]->GetEntries());
            int rangeMin = hCorrSideSideProj[itrigg][iassoc]->GetNbinsX()/2.;
            int rangeMax = hCorrSideSideProj[itrigg][iassoc]->GetNbinsX();
            cout << "side-side : " << hCorrSideSideProj[itrigg][iassoc]->Integral(rangeMin, rangeMax) << "\tfake-fake : " << hCorrBBProj[itrigg][iassoc]->Integral(rangeMin, rangeMax) << endl;
            hCorrSideSideProj[itrigg][iassoc]->Scale(1./hCorrSideSideProj[itrigg][iassoc]->Integral(rangeMin, rangeMax));
            hCorrSideSideProj[itrigg][iassoc]->Draw("HIST E");
            hCorrBBProj[itrigg][iassoc]->Scale(1./hCorrBBProj[itrigg][iassoc]->Integral(rangeMin, rangeMax));
            hCorrBBProj[itrigg][iassoc]->Draw("SAME HIST E");
        }
    }
}
