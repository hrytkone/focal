double Fit(double *x, double *p);
void SetStyle(Bool_t graypalette);
void redrawBorder();

void CheckMissingPionsRatio(TString sInputName = "input.root")
{
    gStyle->SetOptStat(0);
    SetStyle(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    //TH1D *hPionPtFor = (TH1D*)fIn->Get("hPionPtFor");
    //TH1D *hPionPtForDetected = (TH1D*)fIn->Get("hPionPtDetected");
    TH2D *hPionEtaPtRatio = (TH2D*)fIn->Get("hPionEtaPtRatio");

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.58, 0.65, 0.78, 0.85);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.05);

    // ----------------
    // |   Analysis   |
    // ----------------

    /**TCanvas *cCorr = new TCanvas("cCorr", "cCorr", 600, 600);

    TH1D* hRatio = (TH1D*)hPionPtForDetected->Clone("hRatio");
    hRatio->Divide(hPionPtFor);
    hRatio->SetTitle("Ratio of detected #pi^{0}'s to all #pi^{0}'s in FoCal acceptance; p_{T} (GeV); Ratio");
    hRatio->GetYaxis()->SetTitleOffset(1.);
    hRatio->GetYaxis()->SetRangeUser(0., 1.2);

    TF1 *fFit = new TF1("fFit", Fit, 0., 11., 2);
    double par[2];
    hRatio->Fit("fFit","Q0");
    fFit->GetParameters(par);

    std::cout << "Fit function : exp(a/(x + b))" << std::endl;
    std::cout << "\ta = " << par[0] << std::endl;
    std::cout << "\tb = " << par[1] << std::endl;

    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 600, 600);
    cRatio->SetLogx();
    hRatio->Draw();
    fFit->Draw("SAME");**/

    TCanvas *cRatio2D = new TCanvas("cRatio2D", "cRatio2D", 600, 600);
    cRatio2D->cd(0);
    cRatio2D->SetLogy();
    hPionEtaPtRatio->SetTitle("; #eta; p_{T} (GeV)");
    hPionEtaPtRatio->GetYaxis()->SetTitleOffset(0.8);
    hPionEtaPtRatio->GetXaxis()->SetTitleOffset(0.6);
    hPionEtaPtRatio->GetXaxis()->SetTitleSize(0.048);
    hPionEtaPtRatio->GetYaxis()->SetTitleSize(0.048);
    hPionEtaPtRatio->GetXaxis()->SetLabelSize(0.048);
    hPionEtaPtRatio->GetYaxis()->SetLabelSize(0.048);
    hPionEtaPtRatio->GetYaxis()->SetRangeUser(0., 80.);
    hPionEtaPtRatio->Draw("COLZ");
    gStyle->SetPalette(255,0);
    gPad->SetGridy(); gPad->SetGridx();
    gPad->RedrawAxis("f");
    redrawBorder();
}

double Fit(double *x, double *p)
{
    return TMath::Exp(p[0]/(x[0] + p[1]));
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
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadTopMargin(0.05);
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

void redrawBorder()
{
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
   l.SetLineWidth(2);
   l.DrawLine(gPad->GetUxmin(), 80., gPad->GetUxmax(), 80.);
   l.DrawLine(gPad->GetUxmax(), 0.1, gPad->GetUxmax(), 80.);
   l.DrawLine(gPad->GetUxmin(), 0.1, gPad->GetUxmin(), 80.);
   l.DrawLine(gPad->GetUxmin(), 0.1, gPad->GetUxmax(), 0.1);
}
