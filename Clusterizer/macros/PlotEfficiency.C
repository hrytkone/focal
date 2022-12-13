const TString filename = "efficiency_pythiamb_asym-08.root";

const int nPtBin = 6;
double pt[nPtBin+1], limMin = 2, limMax = 20;
double logBW = (log(limMax) - log(limMin))/nPtBin;

const int nEtaBin = 52;
double eta[nEtaBin+1];
double etaBW = 0.05, etamin = 3.2, etamax = 5.8;

TFile *fin;

TH2D *hEtaPtTrue;
TH2D *hPhiEtaTrue;
TH2D *hEtaPtRec;

TH2D *hPhiEta;
TH2D *hPhiTheta;
TH2D *hXY;

TH2D *hPhiEtaGamma;
TH2D *hPhiThetaGamma;
TH2D *hXYGamma;

TH2D *hEtaEff;

TH2D *hEfficiency;
TH1D *hEfficiency_py;
TH1D *hEtaPtTrue_py;
TH1D *hEtaPtRec_py;

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
void SetStyle(Bool_t graypalette);

void PlotEfficiency()
{
    SetStyle(0);

    LoadData();
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    c1->SetLogz();
    hEfficiency->SetTitle("Efficiency ;#eta;p_{T}");
    hEfficiency->Draw("COLZ");

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    c2->SetLogx();
    hEfficiency_py->GetXaxis()->SetRangeUser(2., 20.);
    hEfficiency_py->Rebin(5);
    hEfficiency_py->SetTitle(";p_{T};efficiency");
    hEfficiency_py->Draw("PE");

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
    c3->SetLogx();
    hEtaPtTrue_py->SetTitle("p_{T}^{true} vs p_{T}^{rec};p_{T};efficiency");
    hEtaPtTrue_py->Draw("PE");
    hEtaPtRec_py->Draw("PE SAME");

    TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
    hEPhotonTrue->Draw("PE");
    hEPhotonRec->Draw("PE SAME");

    TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
    hPtMass->RebinY();
    hPtMass->Draw("COLZ");

    // Eta-mass histogram plotting
    TCanvas *c6 = new TCanvas("c6", "c6", 600, 600);

    hEtaMass[0]->Add(hEtaMass[1]);
    hEtaMass[0]->Add(hEtaMass[2]);
    hEtaMass[0]->Add(hEtaMass[3]);
    hEtaMass[0]->Add(hEtaMass[4]);
    hEtaMass[0]->Add(hEtaMass[5]);

    //hEtaMass->RebinY();
    hEtaMass[0]->GetXaxis()->SetRangeUser(3.2, 5.8);
    hEtaMass[0]->GetYaxis()->SetRangeUser(0., 400.);
    hEtaMass[0]->GetXaxis()->SetTitleSize(0.042);
    hEtaMass[0]->GetYaxis()->SetTitleSize(0.042);
    hEtaMass[0]->GetXaxis()->SetLabelSize(0.042);
    hEtaMass[0]->GetYaxis()->SetLabelSize(0.042);
    hEtaMass[0]->SetTitle(";#eta;m_{#gamma#gamma} (MeV/c^{2})");
    //hEtaMass[0]->Fit("pol1");
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
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

    //FitMassDistributions();

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
    hPhiTheta->Draw("same col2 pol");
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
    hPhiThetaGamma->Draw("same col2 pol");
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

    //c10->cd(2);
    //gPad->SetLogz();
    //hXYGamma->SetTitle(";x(cm);y(cm)");
    //hXYGamma->Rebin2D();
    ////hXYGamma->GetXaxis()->SetRangeUser(-49.9,49.9);
    ////hXYGamma->GetYaxis()->SetRangeUser(-49.9,49.9);
    ////hXYGamma->GetXaxis()->SetTitleSize(0.);
    ////hXYGamma->GetYaxis()->SetTitleSize(0.);
    ////hXYGamma->GetXaxis()->SetLabelOffset(999);
    ////hXYGamma->GetXaxis()->SetLabelSize(0.);
    ////hXYGamma->GetYaxis()->SetLabelOffset(999);
    ////hXYGamma->GetYaxis()->SetLabelSize(0.);
    //hXYGamma->Draw("COLZ");
    ////ring21->Draw("same");
    ////ring22->Draw("same");
    ////text2.DrawLatexNDC(0.475, 0.74, "#eta = 3.9");
    ////text2.DrawLatexNDC(0.475, 0.625, "#eta = 4.6");

    c10->SaveAs("xy_markings.png");

    TCanvas *c11 = new TCanvas("c11", "c11", 600, 600);
    hEtaEff->SetTitle(";#eta_{rec};#eta_{true}");
    hEtaEff->Draw("COLZ");
}

void LoadData()
{
    fin = TFile::Open(filename.Data());
    hEtaPtTrue = (TH2D*)fin->Get("hEtaETrue");
    hEtaPtRec = (TH2D*)fin->Get("hEtaERec");
    hEfficiency = (TH2D*)hEtaPtRec->Clone("hEfficiency");
    hEfficiency->Divide(hEtaPtTrue);
    hEtaPtTrue_py = hEtaPtTrue->ProjectionY();
    hEtaPtTrue_py->SetLineColor(kBlack);
    hEtaPtRec_py = hEtaPtRec->ProjectionY();
    hEtaPtRec_py->SetLineColor(kRed);
    hEfficiency_py = (TH1D*)hEtaPtRec_py->Clone("hEfficiency");
    hEfficiency_py->Divide(hEtaPtTrue_py);

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

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetLineScalePS(1);
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
