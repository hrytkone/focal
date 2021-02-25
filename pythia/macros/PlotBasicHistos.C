void PlotBasicHistos(TString sInName = "output.root")
{
    TFile *fIn = TFile::Open(sInName);

    TH1D *hCounter = (TH1D*)fIn->Get("hCounter");
    double nev = hCounter->GetBinContent(1);
    std::cout << "nEvents : " << nev << std::endl;

    TH1D *hPtPi0 = (TH1D*)fIn->Get("hPionPt");
    TH1D *hPtPi0Mid = (TH1D*)fIn->Get("hPionPtMid");
    TH1D *hPtPi0For = (TH1D*)fIn->Get("hPionPtFor");

    hPtPi0->Scale(1.0/nev, "width");
    hPtPi0Mid->Scale(1.0/nev, "width");
    hPtPi0For->Scale(1.0/nev, "width");
    
    gStyle->SetOptStat(0);

    hPtPi0->SetTitle("; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}");
    hPtPi0->SetMarkerStyle(20);
    hPtPi0->SetMarkerSize(.5);
    hPtPi0->SetMarkerColor(1);
    
    hPtPi0Mid->SetMarkerStyle(20);
    hPtPi0Mid->SetMarkerSize(.5);
    hPtPi0Mid->SetMarkerColor(2);
    
    hPtPi0For->SetMarkerStyle(20);
    hPtPi0For->SetMarkerSize(.5);
    hPtPi0For->SetMarkerColor(4);

    TLegend *leg = new TLegend(0.68, 0.68, 0.88, 0.85);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetHeader("PYTHIA #pi^{0} production");
    leg->AddEntry(hPtPi0, "Full #eta", "p");
    leg->AddEntry(hPtPi0Mid, "-1.0 < #eta < 1.0", "p");
    leg->AddEntry(hPtPi0For, "3.2 < #eta < 5.8", "p");
    
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->SetLogy(1);
    hPtPi0->Draw("P");
    hPtPi0Mid->Draw("P SAME");
    hPtPi0For->Draw("P SAME");
    leg->Draw("SAME");

}

