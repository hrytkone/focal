/**
 *  Get parameters for function to emulate the single photon efficiency
 */

double Fit(double *x, double *p);

void PhotonEfficiency()
{
    gStyle->SetOptStat(0);
    
    const int nBins = 75;
    TH1D *hEfficiency = new TH1D("hEfficiency", "hEfficiency", nBins, 0.0, 1500.0);

    for (int i = 1; i <= nBins; i++) hEfficiency->SetBinContent(i, 1.0);
    hEfficiency->SetBinContent(1, 0.695);
    hEfficiency->SetBinContent(2, 0.91);
    hEfficiency->SetBinContent(3, 0.965);
    hEfficiency->SetBinContent(7, 0.99);
    hEfficiency->SetBinContent(25, 0.99);
    hEfficiency->SetBinContent(41, 0.99);

    TF1 *fFit = new TF1("fFit", Fit, 0.0, 1500.0, 1);
    fFit->SetNpx(1000);
    
    double par[1];
    hEfficiency->Fit("fFit", "R+");
    fFit->GetParameters(par);

    std::cout << "Fit result (f(x) = exp(-a/x)) :" << std::endl;
    std::cout << "\ta=" << par[0] << std::endl;

    TCanvas *c = new TCanvas("c", "c", 600, 600);
    hEfficiency->Draw("HIST");
    hEfficiency->GetYaxis()->SetRangeUser(0., 1.1);
    hEfficiency->SetTitle("Single photon efficiency; E (GeV); efficiency");
    fFit->Draw("SAME");
}

double Fit(double *x, double *p)
{
    return TMath::Exp(-p[0]/(x[0]));
}
