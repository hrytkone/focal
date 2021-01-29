void PlotPi0Test(TString sInputName = "output.root")
{

    gStyle->SetOptStat(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hCounter = (TH1D*)fIn->Get("hCounter");
    int nEvent = hCounter->GetBinContent(1);
    std::cout << "Input file contains " << nEvent << " events, proceed to analyse" << std::endl;

    TH2D* hCorrRealPion = (TH2D*)fIn->Get("hCorrRealPions");
    hCorrRealPion->Rebin2D();
    hCorrRealPion->GetYaxis()->SetMaxDigits(3);
    
    TH2D* hCorrFakePion = (TH2D*)fIn->Get("hCorrFakePions");
    hCorrFakePion->Rebin2D();
    hCorrFakePion->GetYaxis()->SetMaxDigits(3);
    
    TH1D* hMassRealPion = (TH1D*)fIn->Get("hMassRealPions");
    hMassRealPion->Rebin();
    hMassRealPion->GetYaxis()->SetMaxDigits(3);
    
    TH1D* hMassFakePion = (TH1D*)fIn->Get("hMassFakePions");
    hMassFakePion->Rebin();
    hMassFakePion->GetYaxis()->SetMaxDigits(3);

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.58, 0.6, 0.78, 0.85);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.05); 

    TLegend *leg2 = new TLegend(0.58, 0.66, 0.78, 0.85);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.05); 

    // ----------------
    // |   Analysis   |
    // ----------------
        
    TCanvas *cMass = new TCanvas("cMass", "cMass", 800, 400);
    cMass->Divide(2, 1);
    
    cMass->cd(1);
    hMassRealPion->SetTitle("Real #pi^{0}; Invariant mass (GeV/c^{2}); Counts");
    hMassRealPion->Draw("HIST");
    
    cMass->cd(2);
    hMassFakePion->SetTitle("Fake #pi^{0}; Invariant mass (GeV/c^{2}); Counts");
    hMassFakePion->Draw("HIST");
    
    TCanvas *cCorr = new TCanvas("cCorr", "cCorr", 800, 400);
    cCorr->Divide(2, 1);
    
    cCorr->cd(1);
    TH1D* hCorrRealPionProj = hCorrRealPion->ProjectionX();
    hCorrRealPionProj->SetTitle("Real #pi^{0}; Invariant mass (GeV/c^{2}); Counts");
    hCorrRealPionProj->Draw("HIST");
    
    cCorr->cd(2);
    TH1D* hCorrFakePionProj = hCorrFakePion->ProjectionX();
    hCorrFakePionProj->SetTitle("All #pi^{0} candidates; Invariant mass (GeV/c^{2}); Counts");
    hCorrFakePionProj->Draw("HIST");

}


// -----------------
// |   Functions   |
// -----------------

