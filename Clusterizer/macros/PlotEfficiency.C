//const TString filename = "efficiency_pi0-gun_MW-120-150.root";
const TString filename = "efficiency_pi0-gun_MW-100-150_matched.root";
//const TString filename = "efficiency_pi0-gun_whole-eta_matched.root";
//const TString filename = "efficiency_pi0-gun_whole-eta_mass-matched.root";

//const TString legHeader = "Mass window [120,150]";
//const TString legHeader = "Mass matched clusters with true p_{T}";
const TString legHeader = "Only clusters matched with mass";

const int nPtBin = 38;
double pt[nPtBin+1], limMin = 2, limMax = 18;
double logBW = (log(limMax) - log(limMin))/nPtBin;

const int nEtaBin = 38;
double eta[nEtaBin+1], etamin = 3.4, etamax = 5.3;
double etaBW = (etamax - etamin)/nEtaBin;

TFile *fin;

TH2D *hEtaPtTrue;
TH2D *hPhiEtaTrue;
TH2D *hEtaPtRec;

TH1D *hEtaPtRec_px;
TH1D *hEtaPtRec_py;
TH1D *hEtaPtTrue_px;
TH1D *hEtaPtTrue_py;

TH2D *hPhiEta;
TH2D *hPhiTheta;
TH2D *hXY;

TH2D *hPhiEtaGamma;
TH2D *hPhiThetaGamma;
TH2D *hXYGamma;

TH2D *hEtaEff;

TH2D *hEfficiency;
TH1D *hEfficiency_py;

TH2D *hPtMass;
TH2D *hEtaMass[nPtBin];

TH1D *hMassCluster[nEtaBin][nPtBin];
TF1 *fFit[nEtaBin][nPtBin];
TF1 *fBg[nEtaBin][nPtBin];
TF1 *fPeak[nEtaBin][nPtBin];

TH1D *hEPhotonTrue;
TH1D *hEPhotonRec;

void LoadData();
void FitMassDistributions();
double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);
void PrintEfficiencyMatrix();
void SetStyle(Bool_t graypalette);

void PlotEfficiency()
{
    SetStyle(0);

    LoadData();
    PrintEfficiencyMatrix();

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    c1->SetLogy();
    hEfficiency->GetZaxis()->SetRangeUser(0., 1.2);
    hEfficiency->SetTitle("Efficiency ;#eta;p_{T}");
    hEfficiency->Draw("COLZ");

    TLegend *legc2 = new TLegend(0.22, 0.2, 0.4, 0.4);
    legc2->SetBorderSize(0); legc2->SetTextSize(0.035);
    legc2->SetHeader(Form("%s", legHeader.Data()));
    legc2->AddEntry(hEtaPtTrue_px, "MC truth", "le");
    legc2->AddEntry(hEtaPtRec_px, "#pi^{0} rec", "le");

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    //hEtaPtTrue_px->GetYaxis()->SetRangeUser(0., 15.);
    hEtaPtTrue_px->Draw("HIST E");
    hEtaPtRec_px->Draw("HIST E SAME");
    legc2->Draw("SAME");

    TLegend *legc3 = new TLegend(0.15, 0.2, 0.5, 0.4);
    legc3->SetBorderSize(0); legc3->SetTextSize(0.035);
    legc3->SetHeader(Form("%s", legHeader.Data()));
    legc3->AddEntry(hEtaPtTrue_px, "MC truth", "le");
    legc3->AddEntry(hEtaPtRec_px, "#pi^{0} rec", "le");

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
    //hEtaPtTrue_py->GetYaxis()->SetRangeUser(3., 5.);
    hEtaPtTrue_py->Draw("HIST E");
    hEtaPtRec_py->Draw("HIST E SAME");
    legc3->Draw("SAME");

    TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
    hPtMass->RebinY();
    hPtMass->GetZaxis()->SetRangeUser(0., 400.);
    hPtMass->Draw("COLZ");

    // Eta-mass histogram plotting
    TCanvas *c6 = new TCanvas("c6", "c6", 600, 600);

    for (int ibin = 1; ibin < nEtaBin; ibin++) {
        hEtaMass[0]->Add(hEtaMass[ibin]);
    }

    //hEtaMass->RebinY();
    //hEtaMass[0]->GetXaxis()->SetRangeUser(4.5, 5.5);
    hEtaMass[0]->GetYaxis()->SetRangeUser(0., 400.);
    hEtaMass[0]->GetXaxis()->SetTitleSize(0.042);
    hEtaMass[0]->GetYaxis()->SetTitleSize(0.042);
    hEtaMass[0]->GetXaxis()->SetLabelSize(0.042);
    hEtaMass[0]->GetYaxis()->SetLabelSize(0.042);
    hEtaMass[0]->SetTitle(";#eta;m_{#gamma#gamma} (MeV/c^{2})");
    //hEtaMass[0]->Fit("pol1");
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    hEtaMass[0]->GetZaxis()->SetRangeUser(0., 400.);
    hEtaMass[0]->Draw("COL");

    TLine *lPionMass = new TLine(3.2, 135., 5.8, 135.);
    lPionMass->SetLineColor(kWhite);
    lPionMass->SetLineWidth(2);
    lPionMass->SetLineStyle(2);
    lPionMass->Draw("SAME");

    //TLine *lEtaMass = new TLine(3.2, 548., 5.8, 548.);
    //lEtaMass->SetLineColor(kWhite);
    //lEtaMass->SetLineWidth(1);
    //lEtaMass->SetLineStyle(9);
    //lEtaMass->Draw("SAME");

    TLatex textPion;
    textPion.SetTextColor(kWhite);
    textPion.SetTextSize(0.045);
    textPion.DrawLatexNDC(.18,.42,"m_{#pi0}=135 MeV/c^{2}");
    c6->SaveAs("pion-mass_updated.pdf");

    //TLatex textEta;
    //textEta.SetTextColor(kWhite);
    //textEta.SetTextSize(0.035);
    //textEta.DrawLatexNDC(.18,.75,"m_{#eta}=548 MeV/c^{2}");

    // Eta-phi histogram plotting
    TCanvas *c7 = new TCanvas("c7", "c7", 1200, 600);
    c7->Divide(2, 1);

    c7->cd(1);
    gPad->SetTheta(90.);
    gPad->SetPhi(0.);
    //gPad->DrawFrame(-5.8, -5.8, 5.8, 5.8);
    //hPhiEtaTrue->RebinY();
    hPhiEtaTrue->GetXaxis()->SetTitleSize(0.042);
    hPhiEtaTrue->GetYaxis()->SetTitleSize(0.042);
    hPhiEtaTrue->GetXaxis()->SetLabelSize(0.042);
    hPhiEtaTrue->GetYaxis()->SetLabelSize(0.042);
    hPhiEtaTrue->SetTitle(";#phi;#eta");
    //hPhiEtaTrue->Draw("lego2 pol");
    //hPhiEtaTrue->Draw("colz psr");
    //hPhiEtaTrue->Draw("same colz pol");
    hPhiEtaTrue->Draw("same colz");


    c7->cd(2);
    gPad->SetTheta(90.);
    gPad->SetPhi(0.);
    //gPad->DrawFrame(-5.8, -5.8, 5.8, 5.8);
    //hPhiEta->RebinY();
    hPhiEta->Divide(hPhiEtaTrue);
    hPhiEta->GetXaxis()->SetTitleSize(0.042);
    hPhiEta->GetYaxis()->SetTitleSize(0.042);
    hPhiEta->GetXaxis()->SetLabelSize(0.042);
    hPhiEta->GetYaxis()->SetLabelSize(0.042);
    hPhiEta->SetTitle(";#phi;#eta");
    //hPhiEta->Draw("lego2 pol");
    //hPhiEta->Draw("same colz pol");
    hPhiEta->Draw("same colz");

    TCanvas *c8 = new TCanvas("c8", "c8", 1200, 600);
    c8->Divide(2, 1);

    c8->cd(1);
    hPhiEta->ProjectionY()->Draw("HIST");

    c8->cd(2);
    hPhiEta->ProjectionX()->Draw("HIST");

    // Theta-phi histogram plotting
    TCanvas *c9 = new TCanvas("c9", "c9", 1200, 600);
    c9->Divide(2,1);

    c9->cd(1);
    gPad->SetTheta(90.);
    gPad->SetPhi(0.);
    //gPad->DrawFrame(-0.049, -0.049, 0.049, 0.049);
    gPad->DrawFrame(-0.07, -0.07, 0.07, 0.07);
    hPhiTheta->SetTitle(";#phi;#theta");
    hPhiTheta->GetYaxis()->SetRangeUser(0.,0.12);
    hPhiTheta->GetXaxis()->SetTitleSize(0.);
    hPhiTheta->GetYaxis()->SetTitleSize(0.);
    hPhiTheta->GetXaxis()->SetLabelOffset(999);
    //hPhiTheta->GetXaxis()->SetLabelSize(0.);
    hPhiTheta->GetYaxis()->SetLabelOffset(999);
    //hPhiTheta->GetYaxis()->SetLabelSize(0.);
    //hPhiTheta->Draw("surf1 pol");
    hPhiTheta->Draw("same col pol");
    //hPhiTheta->Draw("");
    TEllipse *ring1 = new TEllipse(0, 0, 0.04);
    ring1->SetLineColor(kWhite);
    ring1->SetLineWidth(2);
    ring1->SetLineStyle(2);
    ring1->SetFillStyle(0);
    ring1->Draw("same");

    TEllipse *ring2 = new TEllipse(0, 0, 0.02);
    ring2->SetLineColor(kWhite);
    ring2->SetLineWidth(2);
    ring2->SetLineStyle(2);
    ring2->SetFillStyle(0);
    ring2->Draw("same");

    TEllipse *ring3 = new TEllipse(0, 0, 0.0667);
    ring3->SetLineColor(kWhite);
    ring3->SetLineWidth(2);
    ring3->SetLineStyle(2);
    ring3->SetFillStyle(0);
    ring3->Draw("same");

    TLatex text;
    text.SetTextColor(kWhite);
    text.SetTextFont(62);
    text.SetTextSize(0.052);
    text.DrawLatexNDC(0.475, 0.74, "#eta = 3.9");
    text.DrawLatexNDC(0.475, 0.625, "#eta = 4.6");
    gPad->Update();

    c9->cd(2);
    gPad->SetTheta(90.);
    gPad->SetPhi(0.);
    //gPad->DrawFrame(-0.049, -0.049, 0.049, 0.049);
    gPad->DrawFrame(-0.07, -0.07, 0.07, 0.07);
    hPhiThetaGamma->SetTitle(";#phi;#theta");
    hPhiThetaGamma->GetYaxis()->SetRangeUser(0.,0.12);
    hPhiThetaGamma->GetXaxis()->SetTitleSize(0.);
    hPhiThetaGamma->GetYaxis()->SetTitleSize(0.);
    hPhiThetaGamma->GetXaxis()->SetLabelOffset(999);
    //hPhiThetaGamma->GetXaxis()->SetLabelSize(0.);
    hPhiThetaGamma->GetYaxis()->SetLabelOffset(999);
    //hPhiThetaGamma->GetYaxis()->SetLabelSize(0.);
    //hPhiThetaGamma->Draw("surf1 pol");
    hPhiThetaGamma->Draw("same col pol");
    //hPhiThetaGamma->Draw("");

    ring1->Draw("same");
    ring2->Draw("same");
    ring3->Draw("same");
    text.DrawLatexNDC(0.475, 0.74, "#eta = 3.9");
    text.DrawLatexNDC(0.475, 0.625, "#eta = 4.6");
    gPad->Update();

    c9->SaveAs("phi-theta-pol.pdf");

    // XY histogram plotting
    //TCanvas *c10 = new TCanvas("c10", "c10", 1200, 600);
    TCanvas *c10 = new TCanvas("c10", "c10", 600, 600);
    //c10->Divide(2,1);

    c10->cd(1);
    gPad->SetLogz();
    hXY->SetTitle(";x(cm);y(cm)");
    hXY->Rebin2D();
    hXY->GetXaxis()->SetRangeUser(-49.9,49.9);
    hXY->GetYaxis()->SetRangeUser(-49.9,49.9);
    hXY->Draw("COLZ");

    //TEllipse *ring21 = new TEllipse(0, 0, 28.6184);
    TEllipse *ring21 = new TEllipse(0, 0, 42.7151);
    ring21->SetLineColor(kBlack);
    ring21->SetLineWidth(2);
    ring21->SetLineStyle(2);
    ring21->SetFillStyle(0);
    ring21->Draw("same");

    TEllipse *ring22 = new TEllipse(0, 0, 14.2071);
    ring22->SetLineColor(kBlack);
    ring22->SetLineWidth(2);
    ring22->SetLineStyle(2);
    ring22->SetFillStyle(0);
    ring22->Draw("same");

    TLatex text2;
    text2.SetTextColor(kBlack);
    text2.SetTextFont(62);
    text2.SetTextSize(0.052);
    text2.DrawLatexNDC(0.475, 0.74, "#eta = 3.5");
    text2.DrawLatexNDC(0.475, 0.625, "#eta = 4.6");
    gPad->Update();
}

void LoadData()
{
    fin = TFile::Open(filename.Data());
    hEtaPtTrue = (TH2D*)fin->Get("hEtaPtTrue");
    hEtaPtRec = (TH2D*)fin->Get("hEtaPtRec");
    hEtaPtTrue->Scale(1., "width");
    hEtaPtRec->Scale(1., "width");
    //hEfficiency = (TH2D*)hEtaPtTrue->Clone("hEfficiency");
    hEfficiency = (TH2D*)hEtaPtRec->Clone("hEfficiency");
    hEfficiency->Divide(hEtaPtTrue);

    hEtaPtRec_px = hEtaPtRec->ProjectionX("hEtaPtRec_px", 1, nEtaBin-1);
    hEtaPtRec_py = hEtaPtRec->ProjectionY("hEtaPtRec_py", 1, nPtBin-1);
    hEtaPtTrue_px = hEtaPtTrue->ProjectionX("hEtaPtTrue_px", 1, nEtaBin-1);
    hEtaPtTrue_py = hEtaPtTrue->ProjectionY("hEtaPtTrue_py", 1, nPtBin-1);
    hEtaPtRec_px->SetLineColor(kRed);
    hEtaPtRec_py->SetLineColor(kRed);
    hEtaPtTrue_px->SetLineColor(kBlack);
    hEtaPtTrue_py->SetLineColor(kBlack);
    hEtaPtRec_px->GetYaxis()->SetRangeUser(0., 175000.);
    hEtaPtRec_py->GetYaxis()->SetRangeUser(0., 175000.);
    hEtaPtTrue_px->GetYaxis()->SetRangeUser(0., 175000.);
    hEtaPtTrue_py->GetYaxis()->SetRangeUser(0., 175000.);

    hEtaPtRec_px->GetYaxis()->SetTitleSize(0.038);
    hEtaPtRec_py->GetYaxis()->SetTitleSize(0.038);
    hEtaPtTrue_px->GetYaxis()->SetTitleSize(0.038);
    hEtaPtTrue_py->GetYaxis()->SetTitleSize(0.038);
    hEtaPtRec_px->GetXaxis()->SetTitleSize(0.038);
    hEtaPtRec_py->GetXaxis()->SetTitleSize(0.038);
    hEtaPtTrue_px->GetXaxis()->SetTitleSize(0.038);
    hEtaPtTrue_py->GetXaxis()->SetTitleSize(0.038);
    hEtaPtRec_px->GetXaxis()->SetLabelSize(0.038);
    hEtaPtRec_py->GetXaxis()->SetLabelSize(0.038);
    hEtaPtTrue_px->GetXaxis()->SetLabelSize(0.038);
    hEtaPtTrue_py->GetXaxis()->SetLabelSize(0.038);    
    hEtaPtRec_px->GetYaxis()->SetLabelSize(0.038);
    hEtaPtRec_py->GetYaxis()->SetLabelSize(0.038);
    hEtaPtTrue_px->GetYaxis()->SetLabelSize(0.038);
    hEtaPtTrue_py->GetYaxis()->SetLabelSize(0.038);    

    hEtaPtRec_px->SetTitle("; #eta; counts");
    hEtaPtRec_py->SetTitle("; p_{T}; counts");
    hEtaPtTrue_px->SetTitle("; #eta; counts");
    hEtaPtTrue_py->SetTitle("; p_{T}; counts");

    hPtMass = (TH2D*)fin->Get("hPtMass");
    hPhiTheta = (TH2D*)fin->Get("hPhiTheta");
    hPhiEta = (TH2D*)fin->Get("hPhiEta");
    hXY = (TH2D*)fin->Get("hXY");
    hPhiThetaGamma = (TH2D*)fin->Get("hPhiThetaGamma");
    hPhiEtaGamma = (TH2D*)fin->Get("hPhiEtaGamma");
    hXYGamma = (TH2D*)fin->Get("hXYGamma");

    hEtaEff = (TH2D*)fin->Get("hEtaEff");

    hPhiEtaTrue = (TH2D*)fin->Get("hPhiEtaTrue");
    for (int ipt=0; ipt<nPtBin; ipt++)
        hEtaMass[ipt] = (TH2D*)fin->Get(Form("hEtaMass_%d", ipt));

    for (int i=0; i<nEtaBin; i++) {
        for (int j=0; j<nPtBin; j++) {
            hMassCluster[i][j] = (TH1D*)fin->Get(Form("hMassCluster_%d_%d", i, j));
        }
    }

    hEPhotonTrue = (TH1D*)fin->Get("hEPhotonTrue");
    hEPhotonRec = (TH1D*)fin->Get("hEPhotonCluster");
    hEPhotonTrue->SetLineColor(kRed);
    hEPhotonRec->SetLineColor(kBlack);

    cout << "\npt=[ ";
    for (int i = 0; i <= nPtBin; i++) {
        pt[i] = limMin*exp(i*logBW);
        cout << pt[i] << " ";
    }
    cout << "]\neta=[ ";
    for (int i = 0; i <= nEtaBin; i++) {
        eta[i] = etamin + i*etaBW;
        cout << eta[i] << " ";
    }
    cout << "]\n" << endl;
}

void FitMassDistributions()
{
    for (int i=0; i<nEtaBin; i++) {
        for (int j=0; j<nPtBin; j++) {
            fFit[i][j] = new TF1(Form("fFit_%d_%d", i, j), FitFunction, 50, 400, 10);
            fFit[i][j]->SetParameters(0., 0., 0., 0., 0., 1., 1., 135., 5., 10.);
            fFit[i][j]->SetParLimits(7, 130., 160.);
            fFit[i][j]->SetParLimits(8, 0.1, 50.);
            fFit[i][j]->SetParLimits(9, 5., 50.);
            fFit[i][j]->FixParameter(0, 0.);
            fFit[i][j]->FixParameter(1, 0.);
            fFit[i][j]->FixParameter(6, 0.);
            fFit[i][j]->FixParameter(9, 0.);
            fFit[i][j]->SetNpx(1000);

            fPeak[i][j] = new TF1(Form("fPeak_%d_%d", i, j), FitPeak, 50, 400, 5);
            fPeak[i][j]->SetLineColor(kBlue);
            fPeak[i][j]->SetNpx(1000);

            fBg[i][j] = new TF1(Form("fBg_%d_%d", i, j), FitBackground, 50, 400, 5);
            fBg[i][j]->SetLineColor(kBlack);
            fBg[i][j]->SetLineStyle(kDashed);
            fBg[i][j]->SetNpx(1000);

            hMassCluster[i][j]->Fit(Form("fFit_%d_%d", i, j), "SQNR+");

            double chi2 = fFit[i][j]->GetChisquare();
            int ndf = fFit[i][j]->GetNDF();
            cout << "Eta : [ " << eta[i] << " " << eta[i+1] << " ], pT : [ " << pt[j] << " " << pt[j+1] << " ]\t chi2/NDF : " << chi2/ndf << endl;

            //TF1 *fResults = hMassCluster[i][j]->GetFunction(Form("fFit_%d_%d", i, j));
        }
    }
}

double FitPeak(double *x, double *p)
{
    double c1 = p[0];
    double c2 = p[1];
    double mu = p[2];
    double sigma1 = p[3];
    double sigma2 = p[4];

    return c1*TMath::Gaus(x[0], mu, sigma1) + c2*TMath::Gaus(x[0], mu, sigma2);
    //return c1*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma1*sigma1)) + c2*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma2*sigma2));
}

double FitBackground(double *x, double *p)
{
    double a = p[0];
    double b = p[1];
    double c = p[2];
    double d = p[3];
    double e = p[4];

    return a*x[0]*x[0]*x[0]*x[0] + b*x[0]*x[0]*x[0] + c*x[0]*x[0] + d*x[0] + e;
}

double FitFunction(double *x, double *p)
{
    //return FitBackground(x, p) + FitPeak(x, &p[3]);
    return FitBackground(x, p) + FitPeak(x, &p[5]);
}

void PrintEfficiencyMatrix()
{
    cout << "[";
    for (int ieta = 1; ieta <= nEtaBin; ieta++) {
        cout << "[";
        for (int ipt = 1; ipt <= nPtBin; ipt++) {
            TString val = Form("%0.03f", hEfficiency->GetBinContent(ieta, ipt));
            cout << val;
            if (ipt==nPtBin)
                cout << "]," << endl;
            else
                cout << ",";
        }
    }
    cout << "];" << endl;
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetLineScalePS(1);
    //gStyle->SetPalette(kCool);
    if(graypalette) gStyle->SetPalette(8,0);
    //else gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    //gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    //gStyle->SetPadTickX(0);
    //gStyle->SetPadTickY(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(1);
    gStyle->SetFuncColor(kRed);
    gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.042,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.042,"xyz");
    gStyle->SetTitleOffset(0.75,"y");
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleFillColor(kWhite);
    //gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    //gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    //  gStyle->SetFillColor(kWhite);
    gStyle->SetLegendFont(42);
}
