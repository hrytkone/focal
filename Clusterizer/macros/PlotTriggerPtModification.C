#include "include/Filipad.h"
#include "include/rootcommon.h"

const TString filename = "triggPt-mod_binning.root";
//const double eff = 0.329; // pt = [2,3]
const double eff = 0.343; // pt = [3,4]
//const double eff = 0.374; // pt = [4,8]
//const double eff = 0.494; // pt = [8,10]
const int nev = 495000;

const double minpt = 0.;
const double maxpt = 20.;

TFile *fin;
TH1D *hPtTrue;
TH1D *hPtRec;
TH1D *hPtRecCorr;
TH1D *hRatio;

TLegend *leg;

void LoadData();

void PlotTriggerPtModification()
{
    gStyle->SetOptStat(0);
    LoadData();

    TLegend *leg = new TLegend(0.6, 0.53, 0.8, 0.73);
    leg->SetBorderSize(0); leg->SetTextSize(0.038);
    leg->AddEntry(hPtTrue, "MC truth", "lep");
    leg->AddEntry(hPtRec, "Rec.", "lep");
    leg->AddEntry(hPtRecCorr, "Rec. (eff. corrected)", "lep");

    int padID = 0;
    Filipad *fpad = new Filipad(padID, 1.1, 0.35, 100, 100, 0.8, 1);
    fpad->Draw();
    padID++;

    // Upper pad
    TPad *p = fpad->GetPad(1);
    p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
    hset(*hPtTrue, "p_{T} (GeV/c)", "1/N_{ev} dN/dp_{T}", 1.1,1.4, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hPtTrue->Draw("P");
    hPtRec->Draw("P SAME");
    hPtRecCorr->Draw("P SAME");
    leg->Draw("SAME");

    TLatex tl;
    tl.SetTextSize(0.038);
    tl.DrawLatexNDC(0.6, 0.79, "PYTHIA8 simulation");
    tl.DrawLatexNDC(0.6, 0.74, "3.7 < #eta < 5.3");

    // Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    hset( *hRatio, "p_{T} (GeV/c)", "p_{T,rec}/p_{T,true}",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
    hRatio->Draw("P");
}

//************************************************************************************************
//************************************************************************************************

void LoadData()
{
    fin = TFile::Open(filename.Data(), "READ");

    hPtTrue = (TH1D*)fin->Get("hPtTrue");
    hPtRec = (TH1D*)fin->Get("hPtRec");
    hPtRecCorr = (TH1D*)hPtRec->Clone("hPtRecCorr");

    //hPtTrue->Rebin();
    //hPtRec->Rebin();
    //hPtRecCorr->Rebin();

    hPtTrue->SetMarkerColor(kRed);
    hPtTrue->SetMarkerStyle(kFullCircle);
    hPtTrue->SetLineColor(kRed);
    hPtTrue->SetTitle(";p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}");
    hPtTrue->Scale(1./nev);    
    hPtTrue->GetYaxis()->SetTitleOffset(0.8);
    hPtTrue->GetXaxis()->SetRangeUser(minpt,maxpt);
    hPtTrue->GetYaxis()->SetRangeUser(0.00015,0.05);
    hPtTrue->GetYaxis()->SetMaxDigits(2);

    hPtRec->SetMarkerColor(kBlack);
    hPtRec->SetMarkerStyle(kFullSquare);
    hPtRec->SetLineColor(kBlack);

    hPtRecCorr->Scale(1./eff);
    hPtRecCorr->SetMarkerColor(kBlack);
    hPtRecCorr->SetMarkerStyle(kOpenSquare);
    hPtRecCorr->SetLineColor(kBlack);
 
    hPtRec->Scale(1./nev);
    hPtRecCorr->Scale(1./nev);

    hRatio = (TH1D*)hPtRecCorr->Clone("hRatio");
    hRatio->Divide(hPtTrue);
    hRatio->GetXaxis()->SetRangeUser(minpt,maxpt);

    cout << "Integral true : " << hPtTrue->Integral(hPtTrue->FindBin(minpt), hPtTrue->FindBin(maxpt)) << endl;
    cout << "Integral rec. corr. : " << hPtRecCorr->Integral(hPtRecCorr->FindBin(minpt), hPtRecCorr->FindBin(maxpt)) << endl;
    cout << "Ratio : " << hPtTrue->Integral(hPtTrue->FindBin(minpt), hPtTrue->FindBin(maxpt))/hPtRecCorr->Integral(hPtRecCorr->FindBin(minpt), hPtRecCorr->FindBin(maxpt)) << endl;

}