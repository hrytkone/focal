#include "include/Filipad.h"
#include "include/rootcommon.h"

//const TString filename = "efficiency_pythiamb.root";
const TString filename = "efficiency_pi0-gun_asym-08.root";

//const TString legHeader = Form("#splitline{m_{#gamma#gamma} = [100,150] MeV/c^{2}}{3.4 < #eta < 5.3}");
const TString legHeader = Form("#splitline{m_{#gamma#gamma} = [100,150] MeV/c^{2}}{2 < p_{T} < 20 GeV/c}");
//const TString legHeader = "Mass matched clusters with true p_{T}";

//const int nPtBin = 38;
const int nPtBin = 19;
double pt[nPtBin+1], limMin = 2, limMax = 18;
double logBW = (log(limMax) - log(limMin))/nPtBin;

//const int nEtaBin = 38;
const int nEtaBin = 19;
double eta[nEtaBin+1], etamin = 3.4, etamax = 5.3;
double etaBW = (etamax - etamin)/nEtaBin;

TFile *fin;

TH2D *hEtaPtTrue;
TH2D *hPhiEtaTrue;
TH2D *hEtaPtRec;
TH2D *hEtaPtRec_match;
TH2D *hEtaPtRec_match_truept;

TH1D *hEtaPtRec_px;
TH1D *hEtaPtRec_py;
TH1D *hEtaPtRec_match_px;
TH1D *hEtaPtRec_match_py;
TH1D *hEtaPtRec_match_truept_px;
TH1D *hEtaPtRec_match_truept_py;
TH1D *hEtaPtTrue_px;
TH1D *hEtaPtTrue_py;

TH1D *hEtaPtRec_match_px_ratio;
TH1D *hEtaPtRec_match_py_ratio;
TH1D *hEtaPtRec_match_truept_px_ratio;
TH1D *hEtaPtRec_match_truept_py_ratio;

TH2D *hPhiEta;
TH2D *hPhiTheta;
TH2D *hXY;

TH2D *hPhiEtaGamma;
TH2D *hPhiThetaGamma;
TH2D *hXYGamma;

TH2D *hEtaEff;

TH2D *hEfficiency;
TH1D *hEfficiency_py;
TH2D *hEfficiency_match;

TH2D *hPtMass;
TH2D *hEtaMass;

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
void redrawBorder();
void SetStyle(Bool_t graypalette);

void PlotEfficiency()
{
    SetStyle(0);

    LoadData();
    PrintEfficiencyMatrix();

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
    c3->SetLogy();
    hEfficiency->SetContour(50);
    hEfficiency->GetZaxis()->SetRangeUser(0., 1.);
    hEfficiency->SetTitle("Efficiency ;#eta;p_{T} (GeV/c)");
    hEfficiency->Draw("COLZ");

    TLegend *legc2 = new TLegend(0.18, 0.15, 0.4, 0.4);
    legc2->SetBorderSize(0); legc2->SetTextSize(0.035);
    legc2->SetHeader(Form("%s", legHeader.Data()));
    legc2->AddEntry(hEtaPtTrue_px, "MC truth", "le");
    legc2->AddEntry(hEtaPtRec_match_px, "#pi^{0} rec", "le");
    //legc2->AddEntry(hEtaPtRec_match_truept_px, "#pi^{0} rec (true p_{T})", "le");


    Filipad *fpad = new Filipad(0, 1.1, 0.35, 100, 100, 0.8, 1);
    fpad->Draw();

    // Upper pad
    TPad *p = fpad->GetPad(1);
    p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
    hset(*hEtaPtTrue_px, "#eta", "1/N_{ev} dN/d#eta", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hEtaPtTrue_px->GetYaxis()->SetRangeUser(-0.001, 0.03);
    hEtaPtTrue_px->Draw("HIST E");
    hEtaPtRec_match_px->Draw("HIST E SAME");
    //hEtaPtRec_match_truept_px->Draw("HIST E SAME");
    legc2->Draw("SAME");

    // Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    hset( *hEtaPtRec_match_px_ratio, "#eta", "Ratio to MC truth",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
    hEtaPtRec_match_px_ratio->Draw("HIST E");
    //hEtaPtRec_match_truept_px_ratio->Draw("HIST E SAME");

    TLegend *legc3 = new TLegend(0.18, 0.15, 0.4, 0.4);
    legc3->SetBorderSize(0); legc3->SetTextSize(0.035);
    legc3->SetHeader(Form("%s", legHeader.Data()));
    legc3->AddEntry(hEtaPtTrue_py, "MC truth", "le");
    legc3->AddEntry(hEtaPtRec_match_py, "#pi^{0} rec", "le");
    legc3->AddEntry(hEtaPtRec_match_truept_py, "#pi^{0} rec (true p_{T})", "le");

    Filipad *fpad2 = new Filipad(1, 1.1, 0.35, 100, 100, 0.8, 1);
    fpad2->Draw();

    // Upper pad
    TPad *p1 = fpad2->GetPad(1);
    p1->SetTickx(); p1->SetLogx(0); p1->SetLogy(0); p1->cd();
    hset(*hEtaPtTrue_py, "p_{T} (GeV/c)", "1/N_{ev} dN/dp_{T}", 1.1,1.2, 0.05,0.05, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hEtaPtTrue_py->GetYaxis()->SetRangeUser(-0.001, 0.03);
    hEtaPtTrue_py->Draw("HIST E");
    hEtaPtRec_match_py->Draw("HIST E SAME");
    hEtaPtRec_match_truept_py->Draw("HIST E SAME");
    legc3->Draw("SAME");

    // Lower pad
    p1 = fpad2->GetPad(2);
    p1->SetTickx(); p1->SetGridy(1); p1->SetLogx(0), p1->SetLogy(0); p1->cd();
    hset( *hEtaPtRec_match_py_ratio, "p_{T} (GeV/c)", "Ratio to MC truth",1.1,0.7, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
    hEtaPtRec_match_py_ratio->Draw("HIST E");
    hEtaPtRec_match_truept_py_ratio->Draw("HIST E SAME");

    TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
    
    TLine *lPionMassMin1 = new TLine(2., 100., 18., 100.);
    lPionMassMin1->SetLineColor(kBlack);
    lPionMassMin1->SetLineWidth(2);
    lPionMassMin1->SetLineStyle(2);

    TLine *lPionMassMax1 = new TLine(2., 150., 18., 150.);
    lPionMassMax1->SetLineColor(kBlack);
    lPionMassMax1->SetLineWidth(2);
    lPionMassMax1->SetLineStyle(2);
    
    TLatex textPion;
    textPion.SetTextColor(kWhite);
    textPion.SetTextSize(0.045);
    
    hPtMass->GetYaxis()->SetRangeUser(0., 400.);
    hPtMass->GetZaxis()->SetRangeUser(0., 250.);
    hPtMass->GetXaxis()->SetTitleSize(0.042);
    hPtMass->GetYaxis()->SetTitleSize(0.042);
    hPtMass->GetXaxis()->SetLabelSize(0.042);
    hPtMass->GetYaxis()->SetLabelSize(0.042);
    hPtMass->GetZaxis()->SetLabelSize(0.024);
    hPtMass->GetYaxis()->SetTitleOffset(1.4);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.1);    
    hPtMass->SetTitle(";p_{T} (GeV/c);m_{#gamma#gamma} (MeV/c^{2})");
    hPtMass->Draw("COLZ");
    lPionMassMin1->Draw("SAME");
    lPionMassMax1->Draw("SAME");
    redrawBorder();

    // Eta-mass histogram plotting
    TCanvas *c6 = new TCanvas("c6", "c6", 600, 600);

    TLine *lPionMassMin2 = new TLine(etamin, 100., etamax, 100.);
    lPionMassMin2->SetLineColor(kBlack);
    lPionMassMin2->SetLineWidth(2);
    lPionMassMin2->SetLineStyle(2);

    TLine *lPionMassMax2 = new TLine(etamin, 150., etamax, 150.);
    lPionMassMax2->SetLineColor(kBlack);
    lPionMassMax2->SetLineWidth(2);
    lPionMassMax2->SetLineStyle(2);

    hEtaMass->GetYaxis()->SetRangeUser(0., 400.);
    //hEtaMass->GetZaxis()->SetRangeUser(0., 250.);
    hEtaMass->GetXaxis()->SetTitleSize(0.042);
    hEtaMass->GetYaxis()->SetTitleSize(0.042);
    hEtaMass->GetXaxis()->SetLabelSize(0.042);
    hEtaMass->GetYaxis()->SetLabelSize(0.042);
    hEtaMass->GetZaxis()->SetLabelSize(0.024);
    hEtaMass->GetYaxis()->SetTitleOffset(1.4);
    hEtaMass->SetTitle(";#eta;m_{#gamma#gamma} (MeV/c^{2})");
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.1);
    hEtaMass->Draw("COLZ");
    lPionMassMin2->Draw("SAME");
    lPionMassMax2->Draw("SAME");
    redrawBorder();

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

//    TEllipse *ring21 = new TEllipse(0, 0, 42.7151); // eta = 3.5
    TEllipse *ring21 = new TEllipse(0, 0, 24.400); // eta = 3.2
    ring21->SetLineColor(kBlack);
    ring21->SetLineWidth(2);
    ring21->SetLineStyle(2);
    ring21->SetFillStyle(0);
    ring21->Draw("same");

    TEllipse *ring22 = new TEllipse(0, 0, 14.796);
    ring22->SetLineColor(kBlack);
    ring22->SetLineWidth(2);
    ring22->SetLineStyle(2);
    ring22->SetFillStyle(0);
    ring22->Draw("same");

    TEllipse *ring23 = new TEllipse(0, 0, 32.943);
    ring23->SetLineColor(kBlack);
    ring23->SetLineWidth(2);
    ring23->SetLineStyle(2);
    ring23->SetFillStyle(0);
    ring23->Draw("same");

    TLatex text2;
    text2.SetTextColor(kBlack);
    text2.SetTextFont(62);
    text2.SetTextSize(0.052);
    text2.DrawLatexNDC(0.57, 0.78, "#eta = 3.75");
    text2.DrawLatexNDC(0.52, 0.72, "#eta = 4.05");
    text2.DrawLatexNDC(0.47, 0.64, "#eta = 4.55");
    gPad->Update();
}

void LoadData()
{
    fin = TFile::Open(filename.Data());
    hEtaPtTrue = (TH2D*)fin->Get("hEtaPtTrue");
    hEtaPtRec = (TH2D*)fin->Get("hEtaPtRec");
    hEtaPtRec_match = (TH2D*)fin->Get("hEtaPtRec_match");
    hEtaPtRec_match_truept = (TH2D*)fin->Get("hEtaPtRec_match_truept");
    hEfficiency = (TH2D*)hEtaPtRec->Clone("hEfficiency");
    hEfficiency_match = (TH2D*)hEtaPtRec_match->Clone("hEfficiency_match");
    hEfficiency->Divide(hEtaPtTrue);
    hEfficiency_match->Divide(hEtaPtTrue);

    hEtaPtRec_px = hEtaPtRec->ProjectionX("hEtaPtRec_px", 1, nEtaBin);
    hEtaPtRec_py = hEtaPtRec->ProjectionY("hEtaPtRec_py", 1, nPtBin);
    hEtaPtRec_match_px = hEtaPtRec_match->ProjectionX("hEtaPtRec_match_px", 1, nEtaBin);
    hEtaPtRec_match_py = hEtaPtRec_match->ProjectionY("hEtaPtRec_match_py", 1, nPtBin);
    hEtaPtRec_match_truept_px = hEtaPtRec_match_truept->ProjectionX("hEtaPtRec_match_truept_px", 1, nEtaBin);
    hEtaPtRec_match_truept_py = hEtaPtRec_match_truept->ProjectionY("hEtaPtRec_match_truept_py", 1, nPtBin);
    hEtaPtTrue_px = hEtaPtTrue->ProjectionX("hEtaPtTrue_px", 1, nEtaBin);
    hEtaPtTrue_py = hEtaPtTrue->ProjectionY("hEtaPtTrue_py", 1, nPtBin);

    hEtaPtTrue_px->Scale(1., "width");
    hEtaPtTrue_py->Scale(1., "width");
    int c_px = hEtaPtTrue_px->Integral();
    int c_py = hEtaPtTrue_py->Integral();
    hEtaPtTrue_px->Scale(1./c_px);
    hEtaPtTrue_py->Scale(1./c_py);

    hEtaPtRec_px->Scale(1./c_px, "width");
    hEtaPtRec_py->Scale(1./c_py, "width");
    hEtaPtRec_match_px->Scale(1./c_px, "width");
    hEtaPtRec_match_py->Scale(1./c_py, "width");
    hEtaPtRec_match_truept_px->Scale(1./c_px, "width");
    hEtaPtRec_match_truept_py->Scale(1./c_py, "width");

    cout << "Integral px : " << hEtaPtTrue_px->Integral() << endl;
    cout << "Integral py : " << hEtaPtTrue_py->Integral() << endl;

    hEtaPtRec_match_px->SetLineColor(kRed);
    hEtaPtRec_match_py->SetLineColor(kRed);
    hEtaPtRec_match_truept_px->SetLineColor(kBlue);
    hEtaPtRec_match_truept_py->SetLineColor(kBlue);
    hEtaPtTrue_px->SetLineColor(kBlack);
    hEtaPtTrue_py->SetLineColor(kBlack);

    //hEtaPtRec_px->GetYaxis()->SetRangeUser(0., 1.19);
    //hEtaPtRec_py->GetYaxis()->SetRangeUser(0., 1.19);
    //hEtaPtTrue_px->GetYaxis()->SetRangeUser(0., 1.19);
    //hEtaPtTrue_py->GetYaxis()->SetRangeUser(0., 1.19);

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

    hEtaPtRec_px->SetTitle("; #eta; 1/N_{ev} dN/d#eta");
    hEtaPtRec_py->SetTitle("; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}");
    hEtaPtTrue_px->SetTitle("; #eta; 1/N_{ev} dN/d#eta");
    hEtaPtTrue_py->SetTitle("; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}");

    hEtaPtRec_match_px_ratio = (TH1D*)hEtaPtRec_match_px->Clone("hEtaPtRec_match_px_ratio");
    hEtaPtRec_match_py_ratio = (TH1D*)hEtaPtRec_match_py->Clone("hEtaPtRec_match_py_ratio");
    hEtaPtRec_match_truept_px_ratio = (TH1D*)hEtaPtRec_match_truept_px->Clone("hEtaPtRec_match_truept_px_ratio");
    hEtaPtRec_match_truept_py_ratio = (TH1D*)hEtaPtRec_match_truept_py->Clone("hEtaPtRec_match_truept_py_ratio");
    hEtaPtRec_match_px_ratio->Divide(hEtaPtTrue_px);
    hEtaPtRec_match_py_ratio->Divide(hEtaPtTrue_py);
    hEtaPtRec_match_truept_px_ratio->Divide(hEtaPtTrue_px);
    hEtaPtRec_match_truept_py_ratio->Divide(hEtaPtTrue_py);

    hPtMass = (TH2D*)fin->Get("hPtMass_match");
    hEtaMass = (TH2D*)fin->Get("hEtaMass_match");
    //hPtMass = (TH2D*)fin->Get("hPtMass");
    //hEtaMass = (TH2D*)fin->Get("hEtaMass");    

    hPhiTheta = (TH2D*)fin->Get("hPhiTheta");
    hPhiEta = (TH2D*)fin->Get("hPhiEta");
    hXY = (TH2D*)fin->Get("hXY");
    hPhiThetaGamma = (TH2D*)fin->Get("hPhiThetaGamma");
    hPhiEtaGamma = (TH2D*)fin->Get("hPhiEtaGamma");
    hXYGamma = (TH2D*)fin->Get("hXYGamma");
    hEtaEff = (TH2D*)fin->Get("hEtaEff");
    hPhiEtaTrue = (TH2D*)fin->Get("hPhiEtaTrue");

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
            //TString val = Form("%0.03f", hEfficiency->GetBinError(ieta, ipt));
            cout << val;
            if (ipt==nPtBin)
                cout << "]," << endl;
            else
                cout << ",";
        }
    }
    cout << "];" << endl;
}

void redrawBorder()
{
    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.SetLineWidth(2);
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
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
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
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
    gStyle->SetTitleOffset(1.2,"y");
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
