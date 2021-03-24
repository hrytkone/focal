double Fit(double *x, double *p);

void CheckMissingPionsRatio(TString sInputName = "output.root")
{

    gStyle->SetOptStat(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hPionPtFor = (TH1D*)fIn->Get("hPionPtFor");
    TH1D *hPionPtForDetected = (TH1D*)fIn->Get("hPionPtForDetected");

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.58, 0.65, 0.78, 0.85);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.05); 

    // ----------------
    // |   Analysis   |
    // ----------------
    
    TCanvas *cCorr = new TCanvas("cCorr", "cCorr", 600, 600);
    
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
    hRatio->Draw();
    fFit->Draw("SAME");
}

double Fit(double *x, double *p)
{
    return TMath::Exp(p[0]/(x[0] + p[1]));
}
