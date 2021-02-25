double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);

void PlotMassDist(TString sInputName = "output.root")
{

    gStyle->SetOptStat(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------

    TH1D *hMass = (TH1D*)fIn->Get("hMassLeading");
    TH1D *hCounter = (TH1D*)fIn->Get("hCounter");
    int nEvent = hCounter->GetBinContent(1);
    int nTriggPeak = hCounter->GetBinContent(4);
    int nTriggSide = hCounter->GetBinContent(5);
    //hMass->Scale(1./nEvent);

    // ----------------
    // |   Analysis   |
    // ----------------

    TF1 *fFit = new TF1("fFit", FitFunction, 10, 300, 9);
    fFit->SetParameters(0., 0., 0., 0., 0., 135., 135., 30., 60.);
    fFit->SetParLimits(6, 130., 140.);
    //fFit->SetParLimits(8, -10., 10.);
    fFit->SetNpx(1000);

    TF1 *fPeak = new TF1("fPeak", FitPeak, 10, 300, 6);
    fPeak->SetLineColor(kBlue);
    fPeak->SetNpx(1000);

    TF1 *fBg = new TF1("fBg", FitBackground, 10, 300, 3);
    fBg->SetLineColor(kBlack);
    fBg->SetLineStyle(kDashed);
    fBg->SetNpx(1000);

    double par[9];
    TCanvas *cMass = new TCanvas("cMass", "cMass", 600, 600);
    hMass->Fit("fFit", "R");
    fFit->GetParameters(par);
    fBg->SetParameters(par);
    fPeak->SetParameters(&par[3]);
     
    hMass->SetTitle("p_{T} > 3 GeV/c; Invariant mass (MeV/c^{2}); counts");

    hMass->Draw("HIST");
    fFit->Draw("SAME");
    fBg->Draw("SAME");
    fPeak->Draw("SAME");

    // Integrate to get the normalisations
    double sum = fFit->Integral(110, 160);
    double signal = fPeak->Integral(110, 160);
    double bg = fBg->Integral(110, 160);
    double bgside = fBg->Integral(40, 80) + fBg->Integral(210, 280);

    int nTriggFake = (double)(bg/sum)*nTriggPeak;

    std::cout << "\nFit results : " << std::endl;
    std::cout << "\tSignal (peak) : " << signal << std::endl;
    std::cout << "\tBg (peak) : " << bg << std::endl;
    std::cout << "\tBg (side) : " << bgside << std::endl;
    std::cout << "\tS/(S+B) : " << signal/sum << std::endl;
    std::cout << "\tB/(S+B) : " << bg/sum << std::endl;
    std::cout << "\tBg(side)/Bg(peak) : " << bgside/bg << std::endl;
    std::cout << "\nTrigger statistics : " << std::endl;
    std::cout << "\tN trigger (peak) : " << nTriggPeak << std::endl;
    std::cout << "\tN trigger (side) : " << nTriggSide << std::endl;
    std::cout << "\tN trigger (fake) : " << nTriggFake << std::endl;
    std::cout << "\n\tNtriggSide/NtriggFake : " << (double)(bgside/nTriggSide)/(bg/nTriggPeak) << std::endl;

}


// -----------------
// |   Functions   |
// -----------------

double FitPeak(double *x, double *p)
{   
    double c1 = p[0];
    double c2 = p[1];
    double mu1 = p[2];
    double mu2 = p[3];
    double sigma1 = p[4];
    double sigma2 = p[5];

    return c1*TMath::Exp(-(x[0]-mu1)*(x[0]-mu1)/(2*sigma1*sigma1)) + c2*TMath::Exp(-(x[0]-mu2)*(x[0]-mu2)/(2*sigma2*sigma2));
}

double FitBackground(double *x, double *p)
{
    double a = p[0];
    double b = p[1];
    double c = p[2];

    return a*x[0]*x[0] + b*x[0] + c;
}

double FitFunction(double *x, double *p)
{
    return FitBackground(x, p) + FitPeak(x, &p[3]);
}
