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
    
    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 600, 600);
    hRatio->Draw();
}
