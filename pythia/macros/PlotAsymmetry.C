void PlotAsymmetry(TString sInputName = "output.root")
{

    gStyle->SetOptStat(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hEnergyAsymTrue = (TH1D*)fIn->Get("hEnergyAsymTrue");
    TH1D *hEnergyAsymRec  = (TH1D*)fIn->Get("hEnergyAsymRec");
    hEnergyAsymTrue->SetLineColor(kBlack);
    hEnergyAsymTrue->SetLineColor(kRed);

    TCanvas *cAsym = new TCanvas("cAsym", "cAsym", 600, 600);
    hEnergyAsymRec->Draw("HIST E");
    hEnergyAsymTrue->Draw("HIST E SAME");
}
