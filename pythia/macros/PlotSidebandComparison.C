#include "include/Filipad.h"
#include "include/rootcommon.h"

const double xmin = -TMath::Pi()/6;
//const double xmin = -TMath::Pi();
//const double xmin = TMath::Pi()/2.;

const double xmax = TMath::Pi()/3;
//const double xmax = (3./2)*TMath::Pi();

const double alpha = 1;

const int nset = 1;
const int nTriggBins = 2;
const int nAssocBins = 2;
const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

const TString inputname[nset] = {
//    "/home/heimarry/Simulations/focal/analysis_output/2022-12-19_pp-focal_two-sidebands.root",
    //"/home/heimarry/Simulations/focal/analysis_output/2023-01-10_pp-focal_oangle-divided_asym-1.root",
    //"/home/heimarry/Simulations/focal/analysis_output/2023-01-09_pp-focal_oangle-divided_asym-08.root",
    //"/home/heimarry/Simulations/focal/analysis_output/2023-01-10_pp-focal_oangle-divided_asym-05.root"
    "/home/heimarry/Simulations/focal/analysis_output/2023-01-11_pp-focal_thetacut-0002.root"
//    "/home/heimarry/Simulations/focal/analysis_output/2022-12-20_pp-focal_sideband-50-115-160-200.root",
    //"/home/heimarry/Simulations/focal/analysis_output/2022-12-16_pp-focal_sideband-40-80.root",
//    "/home/heimarry/Simulations/focal/analysis_output/2022-12-16_pp-focal_sideband-210-280.root",
//    "/home/heimarry/Simulations/focal/analysis_output/2023-01-09_pp-focal_sideband-160-200.root",
    //"/home/heimarry/Simulations/focal/analysis_output/2022-12-14_pp-focal.root"
};

const TString legHeader = "p-p #sqrt{s} = 14 TeV";
const TString setlabel[nset+1] = {
    "MC truth",
    "[40,80] & [210,280]",
//    "[50,115] & [160,200]",
    //"[40,80]",
//    "[210,280]",
    //"[160,200]"
    //"[300,450]"
};

const EColor cMarker[nset] = {kRed};//, kBlue, kOrange};//, kMagenta};

TFile *fIn[nset];
TH1D *hCounter[nset];
TH2D* hCorrSideSide[nset][nTriggBins][nAssocBins];
TH1D* hCorrSideSideProj[nset][nTriggBins][nAssocBins];
TH2D* hCorrBB[nset][nTriggBins][nAssocBins];
TH1D* hCorrBBProj[nset][nTriggBins][nAssocBins];
TH1D *hSideSideVsBB[nset][nTriggBins][nAssocBins];

TLegend *leg[nTriggBins][nAssocBins];
Filipad *fpad[nTriggBins][nAssocBins];
Filipad *fpadcomp[nTriggBins][nAssocBins];

void LoadData();
void ConfigLegends();
void ConfigHistos();
void DrawFiliPad();

void PlotSidebandComparison()
{
    gStyle->SetOptStat(0);
    LoadData();
    ConfigHistos();
    ConfigLegends();
    DrawFiliPad();
}


// -----------------
// |   Functions   |
// -----------------

void LoadData()
{
    for (int iset = 0; iset < nset; iset++) {
        fIn[iset] = TFile::Open(inputname[iset]);
        hCounter[iset] = (TH1D*)fIn[iset]->Get("hCounter");
        int nEvent = hCounter[iset]->GetBinContent(1);
        std::cout << "Input file " << inputname[iset] << " contains " << nEvent << " events, load histos..." << std::endl;
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;

                hCorrSideSide[iset][itrigg][iassoc] = (TH2D*)fIn[iset]->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                //hCorrSideSide[iset][itrigg][iassoc]->Rebin2D();
                hCorrSideSide[iset][itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);
                hCorrSideSideProj[iset][itrigg][iassoc] = hCorrSideSide[iset][itrigg][iassoc]->ProjectionX();
                hCorrSideSideProj[iset][itrigg][iassoc]->Scale(1./hCorrSideSideProj[iset][itrigg][iassoc]->GetEntries());

                hCorrBB[iset][itrigg][iassoc] = (TH2D*)fIn[iset]->Get(Form("TrueComponents/hCorrBgBg[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                //hCorrBB[iset][itrigg][iassoc]->Rebin2D();
                hCorrBB[iset][itrigg][iassoc]->GetYaxis()->SetMaxDigits(3);
                hCorrBBProj[iset][itrigg][iassoc] = hCorrBB[iset][itrigg][iassoc]->ProjectionX();
                hCorrBBProj[iset][itrigg][iassoc]->Scale(1./hCorrBBProj[iset][itrigg][iassoc]->GetEntries());
            }
        }
    }

    for (int iset = 0; iset < nset; iset++) {
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;
                hSideSideVsBB[iset][itrigg][iassoc] = (TH1D*)hCorrSideSideProj[iset][itrigg][iassoc]->Clone(Form("hSideSideVsBB_%d_[%4.1f,%4.1f][%4.1f,%4.1f]",iset,tlow,tupp,alow,aupp));
                hSideSideVsBB[iset][itrigg][iassoc]->Divide(hCorrBBProj[iset][itrigg][iassoc]);
            }
        }
    }
}

void ConfigLegends()
{
    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

            leg[itrigg][iassoc] = new TLegend(0.55, 0.46, 0.78, 0.71);
            leg[itrigg][iassoc]->SetFillStyle(0); leg[itrigg][iassoc]->SetBorderSize(0); leg[itrigg][iassoc]->SetTextSize(0.05);
            leg[itrigg][iassoc]->SetHeader(Form("#splitline{%s}{[%0.1f,%0.1f][%0.1f,%0.1f]}", legHeader.Data(),triggPt[itrigg],triggPt[itrigg+1],assocPt[iassoc],assocPt[iassoc+1]));
            leg[itrigg][iassoc]->AddEntry(hCorrBBProj[0][itrigg][iassoc], setlabel[0].Data(), "le");
            for (int iset = 0; iset < nset; iset++)
                leg[itrigg][iassoc]->AddEntry(hCorrSideSideProj[iset][itrigg][iassoc], setlabel[iset+1].Data(), "le");
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

                hCorrSideSideProj[iset][itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hCorrSideSideProj[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(xmin, xmax);
                if (cMarker[iset]==kOrange) {
                    hCorrSideSideProj[iset][itrigg][iassoc]->SetLineColor(cMarker[iset]+1);
                    hCorrSideSideProj[iset][itrigg][iassoc]->SetMarkerColor(cMarker[iset]+1);
                } else {
                    hCorrSideSideProj[iset][itrigg][iassoc]->SetLineColor(cMarker[iset]);
                    hCorrSideSideProj[iset][itrigg][iassoc]->SetMarkerColor(cMarker[iset]);
                }
                hCorrSideSideProj[iset][itrigg][iassoc]->SetMarkerStyle(kFullCircle);
                hCorrSideSideProj[iset][itrigg][iassoc]->SetLineWidth(2);
                hCorrSideSideProj[iset][itrigg][iassoc]->SetMarkerSize(0);

                hCorrBBProj[iset][itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hCorrBBProj[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(xmin, xmax);
                hCorrBBProj[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrBBProj[iset][itrigg][iassoc]->SetLineColor(kBlack);
                //hCorrBBProj[iset][itrigg][iassoc]->SetLineStyle(3);
                hCorrBBProj[iset][itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
                hCorrBBProj[iset][itrigg][iassoc]->SetLineWidth(2);
                hCorrBBProj[iset][itrigg][iassoc]->SetMarkerSize(0);

                hSideSideVsBB[iset][itrigg][iassoc]->GetYaxis()->SetTitleOffset(1.);
                hSideSideVsBB[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(xmin, xmax);
                if (cMarker[iset]==kOrange) {
                    hSideSideVsBB[iset][itrigg][iassoc]->SetLineColor(cMarker[iset]+1);
                    hSideSideVsBB[iset][itrigg][iassoc]->SetMarkerColor(cMarker[iset]+1);
                } else {
                    hSideSideVsBB[iset][itrigg][iassoc]->SetLineColor(cMarker[iset]);
                    hSideSideVsBB[iset][itrigg][iassoc]->SetMarkerColor(cMarker[iset]);
                }
                //hSideSideVsBB[iset][itrigg][iassoc]->SetFillColor(cMarker[iset]);
                hSideSideVsBB[iset][itrigg][iassoc]->SetMarkerStyle(kFullCircle);
                hSideSideVsBB[iset][itrigg][iassoc]->SetLineWidth(2);
                hSideSideVsBB[iset][itrigg][iassoc]->SetMarkerSize(0);
            }
        }
    }
}

void DrawFiliPad()
{
    TText *t, *t2;
    t = new TText(.55,0.79,"PYTHIA8 simulation");
    t2 = new TText(.55,0.74,Form("Asymmetry  < %0.1f", alpha));    
    t->SetNDC();
    t2->SetNDC();

    int padID = 0;
    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

            fpadcomp[itrigg][iassoc] = new Filipad(padID, 1.1, 0.35, 100, 100, 0.8, 1);
            fpadcomp[itrigg][iassoc]->Draw();
            padID++;

            // Upper pad
            //int minBin = hCorrSSProj[itrigg][iassoc]->GetMinimumBin();
            //int maxBin = hCorrMassMassProj[itrigg][iassoc]->GetMaximumBin();
            //double rangeMin = 0.5;
            //double rangeMax = hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin) + 10.*hCorrMassMassProj[itrigg][iassoc]->GetBinContent(maxBin);
            TPad *p = fpadcomp[itrigg][iassoc]->GetPad(1);
            p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
            hset(*hCorrBBProj[0][itrigg][iassoc], "#Delta#phi", "1/N dN/d#Delta#phi", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.05,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
            hCorrBBProj[0][itrigg][iassoc]->GetYaxis()->SetRangeUser(0.001, 0.18);
            hCorrBBProj[0][itrigg][iassoc]->Draw("HIST E");
            for (int iset = 0; iset < nset; iset++) hCorrSideSideProj[iset][itrigg][iassoc]->Draw("HIST E SAME");
            //for (int iset = 1; iset < nset; iset++)hCorrBBProj[iset][itrigg][iassoc]->Draw("HIST E SAME");
            hCorrBBProj[0][itrigg][iassoc]->Draw("HIST E SAME");

            leg[itrigg][iassoc]->Draw("SAME");
            t->Draw("SAME");
            t2->Draw("SAME");

            // Lower pad
            p = fpadcomp[itrigg][iassoc]->GetPad(2);
            p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
            hset( *hSideSideVsBB[0][itrigg][iassoc], "#Delta#phi", "f_{side,side}/MC truth",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
            hSideSideVsBB[0][itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.2, 2.2);
            hSideSideVsBB[0][itrigg][iassoc]->Draw("");
            //for (int iset = 0; iset < nset; iset++) {
            //    if (iset==0)
            //        hSideSideVsBB[iset][itrigg][iassoc]->Draw("");
            //    else
            //        hSideSideVsBB[iset][itrigg][iassoc]->Draw("SAME");
            //}
            TLine *l = new TLine(xmin, 1., xmax, 1.);
            l->SetLineStyle(2);
            l->SetLineWidth(2);
            l->Draw("SAME");

            //fpadcomp[itrigg][iassoc]->Print(Form("sb-comp_%d-%d.pdf", itrigg, iassoc));
        }
    }
}
