const int nset = 6;
const int npt = 4;
//const int npt = 10;

//const TString filename = "masses_v5_eta35-55_pt-bins_etacut-02.root";
//const TString filename = "masses_v9.root";
//const TString filename = "masses_pi0-gun_v9.root";
const TString filename = "masses_pythia-2_v9.root";
const TString lowmassfilename = "masses_gamma-gun_v9_2.root";
const TString legEta = "3.4 < #eta < 5.3";
//------------------------------------------------------------------------------

const TString legHeader[npt] = {
    "2 GeV < p_{T} < 3 GeV",
    "3 GeV < p_{T} < 4 GeV",
    "4 GeV < p_{T} < 8 GeV",
    "8 GeV < p_{T} < 16 GeV"
};

//------------------------------------------------------------------------------

const TString legEntry1[nset] = {
    "#alpha < 0.5",
    "#alpha < 0.6",
    "#alpha < 0.7",
    "#alpha < 0.8",
    "#alpha < 0.9",
    "No cut"
};

const double xi[nset] = {0.45, 0.45, 0.45, 0.45, 0.45, 0.45};
const double xf[nset] = {0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
const double yi[nset] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
const double yf[nset] = {0.87, 0.87, 0.87, 0.87, 0.87, 0.87};

const double fitStartingPoint[npt] = {35., 30., 30., 40.};

const double mmin = 50.;
const double mmax = 200.;

//------------------------------------------------------------------------------

TFile *fin;
TFile *finLowMassBg;
TH1D *hCounter;
TH1D *hMassCluster[nset][npt];
TH1D *hMassClusterMixed[nset][npt];
TH1D *hMassClusterBg[nset][npt];
TLegend *leg[nset][npt];
TCanvas *c1[nset];
TCanvas *cStoB;
TCanvas *cEff;
TCanvas *cEffPtFunc;
TLegend *leg2;

TF1 *fFit[nset][npt];
TF1 *fPeak[nset][npt];
TF1 *fBg[nset][npt];

int npion = 0;
double fitPar[10];
double massmin[nset][npt];
double massmax[nset][npt];
double signalToBg[npt][nset] = {0};
double signalToBgErr[npt][nset] = {0};
double peakErr[npt][nset] = {0};
double bgErr[npt][nset] = {0};
double eff[npt][nset] = {0};
double effErr[npt][nset] = {0};

double effPtFunc[nset][npt] = {0};
double effPtFuncErr[nset][npt] = {0};

double peakPos[npt][nset] = {0};
double peakPosErr[npt][nset] = {0};
double asymmetry[nset] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
double asymmteryErr[nset] = {0};
double ptMid[npt] = {2.5, 3.5, 6., 12.};
double ptMidErr[npt] = {0.5, 0.5, 2., 4.};
//double ptMid[npt] = {2.25, 2.75, 3.25, 3.75, 4.5, 5.5, 7, 9, 12.5, 17.5};
//double ptMidErr[npt] = {0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 1., 1., 2.5, 2.5};

const int mColor[npt] = {kAzure, kViolet, kPink, kOrange};//, kSpring, kTeal};//, kBlack, kBlack-1, kBlack-2, kBlack-3};
TGraphErrors *gSignalToBg[npt];
TGraphErrors *gEfficiency[npt];
TGraphErrors *gEffPtFunc[nset];

void LoadData();
void CreateLegends();
void FitMassPeaks();
void CreateGraphs();
double FitPeak(double *x, double *p);
Double_t LorentzianPeak(Double_t *x, Double_t *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);
void SetStyle(Bool_t graypalette);
void redrawBorder();

void PlotMassPeaks()
{
    SetStyle(0);
    //gStyle->SetPalette(kTemperatureMap);
    Int_t palette[6];
    palette[0] = kCyan;
    palette[1] = kCyan+2;
    palette[2] = kMagenta;
    palette[3] = kViolet+1;
    palette[4] = kBlue;
    palette[5] = kMagenta+4;
    gStyle->SetPalette(6, palette);

    LoadData();
    FitMassPeaks();
    CreateGraphs();
    CreateLegends();

    cout << "NPIONS : " << npion << endl;
    for (int j=0; j<nset; j++) {
        c1[j] = new TCanvas(Form("c%d", j), "", 1200, 900);
        c1[j]->Divide(2,2);
        for (int i=0; i<npt; i++) {
            c1[j]->cd(i+1);
            if (i==0 || i==2) {
                gPad->SetRightMargin(0.);
            } else {
                gPad->SetLeftMargin(0.055);
                gPad->SetRightMargin(0.055);
            }

            if (i==0 || i==1) {
                gPad->SetBottomMargin(0.);
                gPad->SetTopMargin(0.06);
            } else {
                gPad->SetTopMargin(0.);
            }

            if (i==2)
                hMassCluster[j][i]->SetTitle(";;");
            else
                hMassCluster[j][i]->SetTitle(";M_{#gamma#gamma}(MeV/c^{2});counts");
            double max = hMassCluster[j][i]->GetBinContent(hMassCluster[j][i]->GetMaximumBin());
            hMassCluster[j][i]->GetXaxis()->SetRangeUser(0.1, 445.);
            hMassCluster[j][i]->GetYaxis()->SetRangeUser(0.1, max+0.14*max);
            //hMassCluster[j][i]->GetYaxis()->SetRangeUser(0., 3000.);
            hMassCluster[j][i]->GetXaxis()->SetTitleSize(0.062);
            hMassCluster[j][i]->GetYaxis()->SetTitleSize(0.062);
            hMassCluster[j][i]->GetXaxis()->SetTitleOffset(0.9);
            hMassCluster[j][i]->GetXaxis()->SetLabelSize(0.052);
            hMassCluster[j][i]->GetYaxis()->SetLabelSize(0.052);
            hMassCluster[j][i]->Draw("HIST E");
            hMassClusterBg[j][i]->Draw("HIST E SAME");
            //hMassClusterMixed[j][i]->Draw("PE SAME");
            fBg[j][i]->Draw("SAME");
            //fPeak[j][i]->Draw("SAME");
            hMassCluster[j][i]->SetLineWidth(2);
            hMassClusterBg[j][i]->SetLineWidth(2);
            hMassClusterBg[j][i]->SetLineColor(kBlue);
            fBg[j][i]->SetLineWidth(2);
            fPeak[j][i]->SetLineWidth(2);
            fFit[j][i]->SetLineWidth(2);
            fFit[j][i]->Draw("SAME");

            TLine *minMassBorder = new TLine(massmin[j][i], 0., massmin[j][i], max);
            TLine *maxMassBorder = new TLine(massmax[j][i], 0., massmax[j][i], max);
            minMassBorder->SetLineStyle(7);
            maxMassBorder->SetLineStyle(7);

            //minMassBorder->Draw("SAME");
            //maxMassBorder->Draw("SAME");
            leg[j][i]->Draw("SAME");
            redrawBorder();
        }
        c1[j]->SaveAs(Form("mass_%d.eps", j));
    }
    cStoB = new TCanvas("cStoB", "", 600, 600);
    cStoB->cd();
    cStoB->SetLogy();
    for (int i=0; i<npt; i++) {
        if (i==0)
            gSignalToBg[i]->Draw("PLA PMC PLC");
        else
            gSignalToBg[i]->Draw("SAME PL PMC PLC");
    }
    leg2->Draw("SAME");
    cStoB->SaveAs("signalToBg.eps");

    cEff = new TCanvas("cEff", "", 600, 600);
    cEff->cd();
    for (int i=0; i<npt; i++) {
        if (i==0)
            gEfficiency[i]->Draw("PLA PMC PLC");
        else
            gEfficiency[i]->Draw("SAME PL PMC PLC");
    }
    leg2->Draw("SAME");
    cEff->SaveAs("efficiency.eps");

    cEffPtFunc = new TCanvas("cEffPtFunc", "", 600, 600);
    cEffPtFunc->cd();
    for (int i=0; i<nset; i++) {
        if (i==0)
            gEffPtFunc[i]->Draw("PLA PMC PLC");
        else
            gEffPtFunc[i]->Draw("SAME PL PMC PLC");
    }
    //leg2->Draw("SAME");
    cEffPtFunc->SaveAs("efficiency-pt.eps");

}

void LoadData()
{
    // To scale the low mass background need to get the areas under both
    // histograms
    double massAreaSidebandLow[nset][npt] = {0};
    double massAreaSidebandHigh[nset][npt] = {0};
    double lowMassArea[nset][npt] = {0};
    double mixedArea[nset][npt] = {0};

    fin = TFile::Open(filename.Data());
    hCounter = (TH1D*)fin->Get("hCounter");
    npion = hCounter->GetEntries();
    for (int iasym=0; iasym<nset; iasym++) {
        for (int ipt=0; ipt<npt; ipt++) {
            hMassCluster[iasym][ipt] = (TH1D*)fin->Get(Form("hMassCluster_%d_%d", iasym, ipt));
            hMassCluster[iasym][ipt]->Rebin(4);
            massAreaSidebandLow[iasym][ipt] = hMassCluster[iasym][ipt]->Integral(hMassCluster[iasym][ipt]->FindBin(0), hMassCluster[iasym][ipt]->FindBin(20));
            massAreaSidebandHigh[iasym][ipt] = hMassCluster[iasym][ipt]->Integral(hMassCluster[iasym][ipt]->FindBin(200), hMassCluster[iasym][ipt]->FindBin(400));
            hMassCluster[iasym][ipt]->Scale(1.);
            hMassCluster[iasym][ipt]->SetMarkerStyle(kDot);
            hMassCluster[iasym][ipt]->SetFillColor(kGray);
            hMassCluster[iasym][ipt]->SetMarkerSize(.5);
            hMassCluster[iasym][ipt]->SetMarkerColor(kBlack);
            hMassCluster[iasym][ipt]->SetLineColor(kBlack);
            hMassCluster[iasym][ipt]->GetYaxis()->SetMaxDigits(3);
            if (ipt==npt-1)
                hMassCluster[iasym][ipt]->GetYaxis()->SetMaxDigits(2);

            hMassClusterMixed[iasym][ipt] = (TH1D*)fin->Get(Form("hMassClusterMixed_%d_%d", iasym, ipt));
            hMassClusterMixed[iasym][ipt]->Rebin(5);
            mixedArea[iasym][ipt] = hMassClusterMixed[iasym][ipt]->Integral(hMassClusterMixed[iasym][ipt]->FindBin(200), hMassClusterMixed[iasym][ipt]->FindBin(400));
            hMassClusterMixed[iasym][ipt]->Scale(massAreaSidebandHigh[iasym][ipt]/mixedArea[iasym][ipt], "width");
            hMassClusterMixed[iasym][ipt]->SetMarkerStyle(kDot);
            hMassClusterMixed[iasym][ipt]->SetMarkerSize(.5);
            hMassClusterMixed[iasym][ipt]->SetMarkerColor(kRed);
            hMassClusterMixed[iasym][ipt]->SetLineColor(kRed);
            hMassClusterMixed[iasym][ipt]->GetYaxis()->SetMaxDigits(3);
        }
    }

    finLowMassBg = TFile::Open(lowmassfilename.Data());
    for (int iasym=0; iasym<nset; iasym++) {
        for (int ipt=0; ipt<npt; ipt++) {
            hMassClusterBg[iasym][ipt] = (TH1D*)finLowMassBg->Get(Form("hMassCluster_%d_%d", iasym, ipt));
            hMassClusterBg[iasym][ipt]->Rebin(4);
            lowMassArea[iasym][ipt] = hMassClusterBg[iasym][ipt]->Integral(hMassClusterBg[iasym][ipt]->FindBin(1), hMassClusterBg[iasym][ipt]->FindBin(20));
            hMassClusterBg[iasym][ipt]->Scale(massAreaSidebandLow[iasym][ipt]/lowMassArea[iasym][ipt]);
            hMassClusterBg[iasym][ipt]->SetMarkerStyle(kDot);
            hMassClusterBg[iasym][ipt]->SetMarkerSize(.5);
            hMassClusterBg[iasym][ipt]->SetMarkerColor(kGray);
            hMassClusterBg[iasym][ipt]->SetLineColor(kBlack);
            //hMassClusterBg[iasym][ipt]->SetFillColor(kGray+1);
            hMassClusterBg[iasym][ipt]->GetYaxis()->SetMaxDigits(3);

            //hMassCluster[iasym][ipt]->Add(hMassClusterBg[iasym][ipt], -1.);
        }
    }
}

void CreateLegends()
{
    for (int iasym=0; iasym<nset; iasym++) {
        for (int ipt=0; ipt<npt; ipt++) {
            leg[iasym][ipt] = new TLegend(xi[ipt], yi[ipt], xf[ipt], yf[ipt]);
            leg[iasym][ipt]->SetFillStyle(0); leg[iasym][ipt]->SetBorderSize(0); leg[iasym][ipt]->SetTextSize(gStyle->GetTextSize()*1.);
            leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], legHeader[ipt].Data(), "");
            leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], legEntry1[iasym], "");
            leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], Form("S/B = %0.03f", signalToBg[ipt][iasym]), "");
            leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], Form("m_{#pi0} = %0.01f#pm%0.01f MeV/c^{2}", peakPos[ipt][iasym], peakPosErr[ipt][iasym]), "");
            leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], "PYTHIA6 + GEANT3", "l");
            leg[iasym][ipt]->AddEntry(hMassClusterBg[iasym][ipt], "#gamma gun + GEANT3", "l");
            //leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], Form("#epsilon = %0.03f#pm%0.03f", eff[ipt][iasym], effErr[ipt][iasym]), "");
            //leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], Form("m_{#pi0} = %0.0f", hMassCluster[iasym][ipt]->GetBinCenter(hMassCluster[iasym][ipt]->GetMaximumBin())), "");
            //leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], Form("N_{entry} = %0.0f", hMassCluster[iasym][ipt]->GetEntries()), "");
            //leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], Form("N_{entry}/N_{#pi0, true} = %0.3f", hMassCluster[iasym][ipt]->GetEntries()/(double)npion), "");
            //leg[iasym][ipt]->AddEntry(hMassCluster[iasym][ipt], legEta.Data(), "");
        }
    }

    // Legend for singal-to-background figure
    leg2 = new TLegend(0.13, 0.7, 0.87, 0.87);
    leg2->SetFillStyle(0); leg2->SetBorderSize(1); leg2->SetTextSize(0.032);
    leg2->SetNColumns(2);
    leg2->SetHeader(legEta.Data(), "C");
    for (int ipt=0; ipt<npt; ipt++)
        leg2->AddEntry(gSignalToBg[ipt], legHeader[ipt].Data(), "pe");
}

void FitMassPeaks()
{
    for (int i=0; i<nset; i++) {
        cout << "\nAsymmetry " << asymmetry[i] << endl;
        for (int j=0; j<npt; j++) {
            fFit[i][j] = new TF1(Form("fFit_%d_%d", i, j), FitFunction, fitStartingPoint[j], 450, 10);
            fFit[i][j]->SetParameters(0., 0., 0., 0., 0., 1., 1., 135., 10., 1.);
            fFit[i][j]->SetParLimits(7, 110., 160.);
            fFit[i][j]->SetParLimits(8, 2., 50.);
            fFit[i][j]->SetParLimits(9, 1., 20.);

            //fFit[i][j] = new TF1(Form("fFit_%d_%d", i, j), FitFunction, fitStartingPoint[j], 400, 8);
            //fFit[i][j]->SetParameters(0., 0., 0., 0., 0., 100., 100., 135.);
            //fFit[i][j]->SetParLimits(8, 120., 160.);

            fFit[i][j]->FixParameter(0, 0.);
            fFit[i][j]->FixParameter(1, 0.);
            //fFit[i][j]->FixParameter(6, 0.);
            //fFit[i][j]->FixParameter(9, 0.);
            fFit[i][j]->SetNpx(1000);

            fPeak[i][j] = new TF1(Form("fPeak_%d_%d", i, j), FitPeak, fitStartingPoint[j], 400, 5);
            fPeak[i][j]->SetLineColor(kBlue);
            fPeak[i][j]->SetNpx(1000);

            fBg[i][j] = new TF1(Form("fBg_%d_%d", i, j), FitBackground, fitStartingPoint[j], 400, 5);
            fBg[i][j]->SetLineColor(kBlack);
            fBg[i][j]->SetLineStyle(kDashed);
            fBg[i][j]->SetNpx(1000);

            TFitResultPtr r = hMassCluster[i][j]->Fit(Form("fFit_%d_%d", i, j), "SQNR+");
            fFit[i][j]->GetParameters(fitPar);
            fBg[i][j]->SetParameters(fitPar);
            fPeak[i][j]->SetParameters(&fitPar[5]);

            //if (fitPar[8] < fitPar[9]) {
            //    massmin[i][j] = fitPar[7]-3.*fitPar[8];
            //    massmax[i][j] = fitPar[7]+3.*fitPar[8];
            //} else {
            //    massmin[i][j] = fitPar[7]-6.*fitPar[9];
            //    massmax[i][j] = fitPar[7]+6.*fitPar[9];
            //}

            massmin[i][j] = mmin;
            massmax[i][j] = mmax;

            peakPos[j][i] = fitPar[7];
            peakPosErr[j][i] = fFit[i][j]->GetParError(7);

            if (!r)
                bgErr[j][i] = fBg[i][j]->IntegralError(massmin[i][j], massmax[i][j], r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());
            else
                bgErr[j][i] = 0.;

            //cout << "\nHistogram : \t" << hMassCluster[i][j]->Integral(hMassCluster[i][j]->FindBin(110), hMassCluster[i][j]->FindBin(160)) << endl;
            //cout << "Fit : \t" << fFit[i][j]->Integral(110, 160)/hMassCluster[i][j]->GetBinWidth(0) << endl;
            //cout << "Fit (peak) : \t" << fPeak[i][j]->Integral(110, 160)/hMassCluster[i][j]->GetBinWidth(0)<< endl;
            //cout << "Fit (bg) : \t" << fBg[i][j]->Integral(110, 160)/hMassCluster[i][j]->GetBinWidth(0)<< endl;
            //cout << "true : \t" << hCounter->GetBinContent(j+1) << endl;

            //cout << legHeader[j] << "  Nrec=" << fPeak[i][j]->Integral(110, 160) << "\tNtrue=" << hCounter->GetBinContent(j+1) << endl;
            //cout << legHeader[j] << "  Nrec=" << hMassCluster[i][j]->Integral(hMassCluster[i][j]->FindBin(110), hMassCluster[i][j]->FindBin(160)) << "\tNtrue=" << hCounter->GetBinContent(j+1) << "\tratio=" <<  hMassCluster[i][j]->Integral(hMassCluster[i][j]->FindBin(110), hMassCluster[i][j]->FindBin(160))/hCounter->GetBinContent(j+1) << endl;
            cout << legHeader[j] << "  Nrec=" << hMassCluster[i][j]->Integral(hMassCluster[i][j]->FindBin(massmin[i][j]), hMassCluster[i][j]->FindBin(massmax[i][j])) - fBg[i][j]->Integral(massmin[i][j], massmax[i][j])/hMassCluster[i][j]->GetBinWidth(0) << "\tNtrue=" << hCounter->GetBinContent(j+1) << "\tsigma=" << fitPar[8] << "\tpeak=" << fitPar[7] << endl;
            //cout << legHeader[j] << "  Nrec=" << hMassCluster[i][j]->GetEntries() << "\tNtrue=" << hCounter->GetBinContent(j+1) << "\tsigma=" << fitPar[8] << "\tpeak=" << fitPar[7] << endl;
        }
    }
}

void CreateGraphs()
{
    for (int i=0; i<npt; i++) {
        for (int j=0; j<nset; j++) {
            //signalToBg[i][j] = fPeak[j][i]->Integral(massmin[i][j], massmax[i][j])/fFit[j][i]->Integral(massmin[i][j], massmax[i][j]);
            //signalToBg[i][j] = fPeak[j][i]->Integral(110, 160)/fBg[j][i]->Integral(110, 160);
            //signalToBg[i][j] = (hMassCluster[j][i]->Integral(hMassCluster[j][i]->FindBin(massmin[i][j]), hMassCluster[j][i]->FindBin(massmax[i][j])) - fBg[j][i]->Integral(massmin[i][j], massmax[i][j])/hMassCluster[j][i]->GetBinWidth(0))/(fFit[j][i]->Integral(massmin[i][j], massmax[i][j])/hMassCluster[j][i]->GetBinWidth(0));
            double intAll = hMassCluster[j][i]->IntegralAndError(hMassCluster[j][i]->FindBin(massmin[j][i]), hMassCluster[j][i]->FindBin(massmax[j][i]), peakErr[j][i]);
            double intBg = fBg[j][i]->Integral(massmin[j][i], massmax[j][i])/hMassCluster[j][i]->GetBinWidth(0);
            signalToBg[i][j] = intAll/intBg - 1.;
            cout << bgErr[i][j] << endl;
            signalToBgErr[i][j] = signalToBg[i][j]*TMath::Sqrt(peakErr[i][j]*peakErr[i][j]/(intAll*intAll) + bgErr[i][j]*bgErr[i][j]/(intBg*intBg));
        }
        gSignalToBg[i] = new TGraphErrors(nset, asymmetry, signalToBg[i], 0, signalToBgErr[i]);
        gSignalToBg[i]->SetTitle(";#alpha;S/B");
        gSignalToBg[i]->SetMarkerColor(mColor[i]);
        gSignalToBg[i]->SetLineColor(mColor[i]);
        gSignalToBg[i]->GetYaxis()->SetRangeUser(0.05, 100.);
        gSignalToBg[i]->GetYaxis()->SetTitleSize(0.052);
        gSignalToBg[i]->GetXaxis()->SetTitleSize(0.052);
        gSignalToBg[i]->GetYaxis()->SetLabelSize(0.042);
        gSignalToBg[i]->GetXaxis()->SetLabelSize(0.042);
        gSignalToBg[i]->GetYaxis()->SetTitleOffset(0.9);
        gSignalToBg[i]->GetXaxis()->SetTitleOffset(0.8);
        gSignalToBg[i]->SetMarkerStyle(20);
        //gSignalToBg[i]->SetMarkerSize(0.5);

        cout << "\npt = " << legHeader[i]<< endl;
        for (int j=0; j<nset; j++) {
            //eff[i][j] = fPeak[j][i]->Integral(110, 160)/hCounter->GetBinContent(i+1);
            if (i==npt-1) {
                eff[i][j] = hMassCluster[j][i]->Integral(hMassCluster[j][i]->FindBin(massmin[j][i]), hMassCluster[j][i]->FindBin(massmax[j][i]))/(hCounter->GetBinContent(i+1));
                effPtFunc[j][i] = hMassCluster[j][i]->Integral(hMassCluster[j][i]->FindBin(massmin[j][i]), hMassCluster[j][i]->FindBin(massmax[j][i]))/(hCounter->GetBinContent(i+1));
            } else {
                eff[i][j] = (hMassCluster[j][i]->Integral(hMassCluster[j][i]->FindBin(massmin[j][i]), hMassCluster[j][i]->FindBin(massmax[j][i])) - fBg[j][i]->Integral(massmin[j][i], massmax[j][i])/hMassCluster[j][i]->GetBinWidth(0))/(hCounter->GetBinContent(i+1));
                effPtFunc[j][i] = (hMassCluster[j][i]->Integral(hMassCluster[j][i]->FindBin(massmin[j][i]), hMassCluster[j][i]->FindBin(massmax[j][i])) - fBg[j][i]->Integral(massmin[j][i], massmax[j][i])/hMassCluster[j][i]->GetBinWidth(0))/(hCounter->GetBinContent(i+1));
            }
            cout << "asymmetry " << asymmetry[j] << "\tEff : " << eff[i][j] << "\tS/B : " << signalToBg[i][j] << "\terr : " << signalToBgErr[i][j] << endl;

            //eff[i][j] = hMassCluster[j][i]->Integral(hMassCluster[j][i]->FindBin(110), hMassCluster[j][i]->FindBin(160))/hCounter->GetBinContent(i+1);
            //double dPeak = fFit[j][i]->IntegralError(110, 160)/fFit[j][i]->Integral(110, 160);
            //double dTrue = hCounter->GetBinError(i+1)/hCounter->GetBinContent(i+1);
            effErr[i][j] = TMath::Sqrt(peakErr[j][i]*peakErr[j][i] + bgErr[j][i]*bgErr[j][i])/hCounter->GetBinContent(i+1);
            effPtFuncErr[j][i] = TMath::Sqrt(peakErr[j][i]*peakErr[j][i] + bgErr[j][i]*bgErr[j][i])/hCounter->GetBinContent(i+1);
        }
        gEfficiency[i] = new TGraphErrors(nset, asymmetry, eff[i], asymmteryErr, effErr[i]);
        gEfficiency[i]->SetTitle(";#alpha;N_{#pi0}^{rec}/N_{#pi0}^{true}");
        gEfficiency[i]->SetMarkerColor(mColor[i]);
        gEfficiency[i]->SetLineColor(mColor[i]);
        gEfficiency[i]->GetYaxis()->SetRangeUser(0., 1.5);
        gEfficiency[i]->GetYaxis()->SetTitleSize(0.052);
        gEfficiency[i]->GetXaxis()->SetTitleSize(0.052);
        gEfficiency[i]->GetYaxis()->SetLabelSize(0.042);
        gEfficiency[i]->GetXaxis()->SetLabelSize(0.042);
        gEfficiency[i]->GetYaxis()->SetTitleOffset(1.1);
        gEfficiency[i]->GetXaxis()->SetTitleOffset(0.8);
        gEfficiency[i]->SetMarkerStyle(20);
        //gEfficiency[i]->SetMarkerSize(0.5);
    }

    for (int i=0; i<nset; i++) {
        gEffPtFunc[i] = new TGraphErrors(npt, ptMid, effPtFunc[i], ptMidErr, effPtFuncErr[i]);
        gEffPtFunc[i]->SetTitle(";p_{T}(GeV/c);N_{#pi0}^{rec}/N_{#pi0}^{true}");
        gEffPtFunc[i]->SetMarkerColor(mColor[i]);
        gEffPtFunc[i]->SetLineColor(mColor[i]);
        gEffPtFunc[i]->GetYaxis()->SetRangeUser(0., 1.);
        gEffPtFunc[i]->SetMarkerStyle(20);
        //gEffPtFunc[i]->SetMarkerSize(0.5);
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

// Lorentzian Peak function
Double_t LorentzianPeak(Double_t *x, Double_t *p) {
    double c = p[0];
    double gamma = p[1];
    double mu = p[2];

    return (0.5*c*gamma/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-mu)*(x[0]-mu) + .25*gamma*gamma);
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
    return FitBackground(x, p) + FitPeak(x, &p[5]);
    //return FitBackground(x, p) + LorentzianPeak(x, &p[5]);
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineScalePS(1);
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
    gStyle->SetHistLineWidth(2);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kRed);
    gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.042,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    //gStyle->SetTitleSize(0.042,"xyz");
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

void redrawBorder()
{
    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.SetLineWidth(2);
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymin());
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmin(), gPad->GetUymax());
    gPad->Update();
}