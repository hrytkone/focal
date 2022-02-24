double Fit(double *x, double *p);
void SetStyle(Bool_t graypalette);

void PlotEfficiency(TString sInputName = "input.root")
{
    gStyle->SetOptStat(0);
    SetStyle(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hPionPtTrue = (TH1D*)fIn->Get("hPionPt");
    hPionPtTrue->Rebin();
    TH1D *hPionPtRec = (TH1D*)fIn->Get("hRecPionPt");
    hPionPtRec->Rebin();

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.3, 0.25, 0.6, 0.42);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.035);

    // ----------------
    // |   Analysis   |
    // ----------------
    TH1D* hRatio = (TH1D*)hPionPtRec->Clone("hRatio");
    hRatio->Divide(hPionPtTrue);
    hRatio->SetTitle("; p_{T} (GeV); Efficiency");
    hRatio->GetYaxis()->SetTitleOffset(1.);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerSize(0.5);
    hRatio->GetYaxis()->SetRangeUser(0., 1.19);
    hRatio->GetXaxis()->SetRangeUser(0.1, 12);

    TF1 *fFit = new TF1("fFit", Fit, 2., 12., 1);
    fFit->SetParameter(0, 0.98);
    double par, parErr;
    hRatio->Fit("fFit","RQ0");
    par = fFit->GetParameter(0);
    parErr = fFit->GetParError(0);

    std::cout << "Fit function : constant" << std::endl;
    std::cout << "\ta = " << par << " +- " << parErr << std::endl;

    leg1->AddEntry(hRatio, "PYTHIA8 simulation", "ep");
    leg1->AddEntry(fFit, Form("Fit constant : %0.04f #pm %0.04f", par, parErr), "l");

    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 600, 600);
    hRatio->Draw("P");
    fFit->Draw("SAME");
    leg1->Draw("SAME");

}

double Fit(double *x, double *p)
{
    //return TMath::Exp(p[0]/(x[0] + p[1]));
    //return p[0]*x[0] + p[1];
    return p[0];
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    //gStyle->Reset("Plain");
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineScalePS(1);
    //if(graypalette) gStyle->SetPalette(8,0);
    //else gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    //gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    //gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.035,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    //gStyle->SetTitleSize(0.035,"xyz");
    //gStyle->SetTitleOffset(1.25,"y");
    //gStyle->SetTitleOffset(1.2,"x");
    //gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    //  gStyle->SetFillColor(kWhite);
    gStyle->SetLegendFont(42);
}
