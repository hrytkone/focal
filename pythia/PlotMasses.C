const int nPi0PtBins = 15;
double pi0Pt[nPi0PtBins+1] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 50.0};

void PlotMasses(TString sInputName = "output.root")
{
    TFile *fIn = TFile::Open(sInputName);

    TH1D *hMassMid[nPi0PtBins];
    TH1D *hMassFor[nPi0PtBins];

    for (int i=0; i<nPi0PtBins; i++) {
        hMassMid[i] = (TH1D*)fIn->Get(Form("hPi0MassMid%d", i));
        hMassFor[i] = (TH1D*)fIn->Get(Form("hPi0MassFor%d", i));
    }

    gStyle->SetOptStat(0);

    hMassMid[0]->SetLineColor(2);
    hMassFor[0]->SetLineColor(4);
    
    TLegend *leg = new TLegend(0.48, 0.68, 0.78, 0.85);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
    leg->AddEntry(hMassMid[0], "-1.0 < #eta < 1.0", "l");
    leg->AddEntry(hMassFor[0], "3.2 < #eta < 5.8", "l");
    
    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 900);
    c1->Divide(5, 3);

    for (int i=0; i<nPi0PtBins; i++) {
        c1->cd(i+1);
        
        hMassFor[i]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c; m_{#gamma#gamma} (MeV/c^{2}); ", pi0Pt[i], pi0Pt[i+1]));
        hMassFor[i]->SetMarkerStyle(20);
        hMassFor[i]->SetLineColor(2);
        //hCorrMidProj->GetYaxis()->SetRangeUser(0.01, 2.0);
    
        hMassMid[i]->SetMarkerStyle(20);
        hMassMid[i]->SetLineColor(4);
        
        hMassFor[i]->Draw("HIST");
        hMassMid[i]->Draw("HIST SAME");
        if (i==0) leg->Draw("SAME");
    }
}
