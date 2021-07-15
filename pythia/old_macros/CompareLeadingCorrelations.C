TString sInputNames[2] = {"2021-02-03_test-data/1400m_pThard-1_corr-ntrigg.root", "2021-02-03_test-data/1000m_pThard-1_leading-ntrigg.root"};

void CompareLeadingCorrelations()
{

    gStyle->SetOptStat(0);

    TFile *fIn1 = TFile::Open(sInputNames[0]);
    TFile *fIn2 = TFile::Open(sInputNames[1]);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hCounter1 = (TH1D*)fIn1->Get("hCounter");
    TH1D *hCounter2 = (TH1D*)fIn2->Get("hCounter");
    
    int nEvent1 = hCounter1->GetBinContent(1);
    int nEvent2 = hCounter2->GetBinContent(1);
    std::cout << "Input file 1 contains " << nEvent1 << " events, input file 2 " << nEvent2 << " events" << std::endl;
    
    int nTriggReal = hCounter1->GetBinContent(4);
    int nTriggFake = hCounter1->GetBinContent(5);
    int nTriggPeak = 0.807*hCounter2->GetBinContent(4);
    int nTriggSide = hCounter2->GetBinContent(5);
    std::cout << "\nNumber of triggers: " << std::endl;
    std::cout << "\treal : " << nTriggReal << "\tper event : " << (double)nTriggReal/nEvent1 << std::endl;
    std::cout << "\tfake : " << nTriggFake << "\tper event : " << (double)nTriggFake/nEvent1 << std::endl;
    std::cout << "\tpeak : " << nTriggPeak << "\tper event : " << (double)nTriggPeak/nEvent2 << std::endl;
    std::cout << "\tside : " << nTriggSide << "\tper event : " << (double)nTriggSide/nEvent2 << std::endl;

    TH2D *hCorrSS = (TH2D*)fIn1->Get("hCorrSS"); hCorrSS->Rebin(4);
    TH2D *hCorrSB = (TH2D*)fIn1->Get("hCorrSB"); hCorrSB->Rebin(4);
    TH2D *hCorrBS = (TH2D*)fIn1->Get("hCorrBS"); hCorrBS->Rebin(4);
    TH2D *hCorrBB = (TH2D*)fIn1->Get("hCorrBB"); hCorrBB->Rebin(4);
    
    TH2D *hCorrReal = (TH2D*)fIn2->Get("hCorrFor0:0"); hCorrReal->Rebin(4);
    TH2D *hCorrMassMass = (TH2D*)fIn2->Get("hCorrMassMass0:0"); hCorrMassMass->Rebin(4);
    TH2D *hCorrMassSide = (TH2D*)fIn2->Get("hCorrMassSide0:0"); hCorrMassSide->Rebin(4);
    TH2D *hCorrSideMass = (TH2D*)fIn2->Get("hCorrSideMass0:0"); hCorrSideMass->Rebin(4);
    TH2D *hCorrSideSide = (TH2D*)fIn2->Get("hCorrSideSide0:0"); hCorrSideSide->Rebin(4);

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.58, 0.7, 0.78, 0.85);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.035); 
    
    TLegend *leg2 = new TLegend(0.58, 0.7, 0.78, 0.85);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.035); 
    
    TLegend *leg3 = new TLegend(0.58, 0.7, 0.78, 0.85);
    leg3->SetFillStyle(0); leg3->SetBorderSize(0); leg3->SetTextSize(0.035); 
    
    TLegend *leg4 = new TLegend(0.58, 0.7, 0.78, 0.85);
    leg4->SetFillStyle(0); leg4->SetBorderSize(0); leg4->SetTextSize(0.035); 

    // ----------------
    // |   Analysis   |
    // ----------------

    // "measured" correlations
    TH1D *hCorrMassMassProj = hCorrMassMass->ProjectionX();
    hCorrMassMassProj->SetLineColor(kBlack);
    hCorrMassMassProj->SetTitle("f_{mass,mass} comparison");
    hCorrMassMassProj->Scale(1./nEvent2);
    hCorrMassMassProj->Scale(1./nTriggPeak);
            
    TH1D *hCorrMassSideProj = hCorrMassSide->ProjectionX();
    hCorrMassSideProj->SetLineColor(kBlack);
    hCorrMassSideProj->SetTitle("f_{mass,side} comparison");
    hCorrMassSideProj->Scale(1./nEvent2);
    hCorrMassSideProj->Scale(1./nTriggPeak);
            
    TH1D *hCorrSideMassProj = hCorrSideMass->ProjectionX();
    hCorrSideMassProj->SetLineColor(kBlack);
    hCorrSideMassProj->SetTitle("f_{side,mass} comparison");
    hCorrSideMassProj->Scale(1./nEvent2);
    hCorrSideMassProj->Scale(1./nTriggSide);
            
    TH1D *hCorrSideSideProj = hCorrSideSide->ProjectionX();
    hCorrSideSideProj->SetLineColor(kBlack);
    hCorrSideSideProj->SetTitle("f_{side,side} comparison");
    hCorrSideSideProj->Scale(1./nEvent2);
    hCorrSideSideProj->Scale(1./nTriggSide);

    // real correlations
    TH1D *hCorrSSProj = hCorrSS->ProjectionX();
    hCorrSSProj->Scale(1./nEvent1);
    hCorrSSProj->Scale(1./nTriggReal);
    
    TH1D *hCorrSBProj = hCorrSB->ProjectionX();
    hCorrSBProj->Scale(1./nEvent1);
    hCorrSBProj->Scale(1./nTriggReal);
    
    TH1D *hCorrBSProj = hCorrBS->ProjectionX();
    hCorrBSProj->Scale(1./nEvent1);
    hCorrBSProj->Scale(1./nTriggFake);
    
    TH1D *hCorrBBProj = hCorrBB->ProjectionX();
    hCorrBBProj->Scale(1./nEvent1);
    hCorrBBProj->Scale(1./nTriggFake);
    hCorrBBProj->SetLineColor(kRed);
    
    // "measured" correlations by summing real correlations
    TH1D *hCorrMassMassRec = (TH1D*)hCorrSSProj->Clone("hCorrSSrec");
    hCorrMassMassRec->Add(hCorrSBProj);
    hCorrMassMassRec->Add(hCorrBSProj);
    hCorrMassMassRec->Add(hCorrBBProj);
    hCorrMassMassRec->SetLineColor(kRed);
    
    TH1D *hCorrMassSideRec = (TH1D*)hCorrSBProj->Clone("hCorrSBrec");
    hCorrMassSideRec->Add(hCorrBBProj);
    hCorrMassSideRec->SetLineColor(kRed);
    
    TH1D *hCorrSideMassRec = (TH1D*)hCorrBSProj->Clone("hCorrBSrec");
    hCorrSideMassRec->Add(hCorrBBProj);
    hCorrSideMassRec->SetLineColor(kRed);
    
    // Ratios
    TH1D *hRatioMassMass = (TH1D*)hCorrMassMassProj->Clone("hRatioMassMass");
    hRatioMassMass->Divide(hCorrMassMassRec);
    hRatioMassMass->SetMarkerStyle(7);
    hRatioMassMass->SetTitle("f_{mass,mass}/(f_{SS} + f_{SB} + f_{BS} + f_{BB})");
    
    TH1D *hRatioMassSide = (TH1D*)hCorrMassSideProj->Clone("hRatioMassSide");
    hRatioMassSide->Divide(hCorrMassSideRec);
    hRatioMassSide->SetMarkerStyle(7);
    hRatioMassSide->SetTitle("f_{mass,side}/(f_{SB} + f_{BB})");
    
    TH1D *hRatioSideMass = (TH1D*)hCorrSideMassProj->Clone("hRatioSideMass");
    hRatioSideMass->Divide(hCorrSideMassRec);
    hRatioSideMass->SetMarkerStyle(7);
    hRatioSideMass->SetTitle("f_{side,mass}/(f_{BS} + f_{BB})");
    
    TH1D *hRatioSideSide = (TH1D*)hCorrSideSideProj->Clone("hRatioSideSide");
    hRatioSideSide->Divide(hCorrBBProj);
    hRatioSideSide->SetMarkerStyle(7);
    hRatioSideSide->SetTitle("f_{side,side}/f_{BB}");

    leg1->AddEntry(hCorrMassMassProj, "f_{mass,mass}", "l");
    leg1->AddEntry(hCorrMassMassRec, "f_{SS} + f_{SB} + f_{BS} + f_{BB}", "l");
    leg2->AddEntry(hCorrMassSideProj, "f_{mass,side}", "l");
    leg2->AddEntry(hCorrMassSideRec, "f_{SB} + f_{BB}", "l");
    leg3->AddEntry(hCorrSideMassProj, "f_{side,mass}", "l");
    leg3->AddEntry(hCorrSideMassRec, "f_{BS} + f_{BB}", "l");
    leg4->AddEntry(hCorrSideSideProj, "f_{side,side}", "l");
    leg4->AddEntry(hCorrBBProj, "f_{BB}", "l");

    
    // Draw histograms
    TCanvas *cCorrMassMassComp = new TCanvas("cCorrMassMassComp", "cCorrMassMassComp", 600, 600);
    cCorrMassMassComp->SetLogy();
    hCorrMassMassProj->Draw("HIST");
    hCorrMassMassRec->Draw("HIST SAME");
    //hCorrMassMassProj->GetYaxis()->SetRangeUser(50, 150000.0);
    leg1->Draw("SAME");
    
    TCanvas *cCorrMassSideComp = new TCanvas("cCorrMassSideComp", "cCorrMassSideComp", 600, 600);
    cCorrMassSideComp->SetLogy();
    hCorrMassSideProj->Draw("HIST");
    hCorrMassSideRec->Draw("HIST SAME");
    //hCorrMassSideProj->GetYaxis()->SetRangeUser(10, 150000.0);
    leg2->Draw("SAME");
    
    TCanvas *cCorrSideMassComp = new TCanvas("cCorrSideMassComp", "cCorrSideMassComp", 600, 600);
    cCorrSideMassComp->SetLogy();
    hCorrSideMassProj->Draw("HIST");
    //hCorrSideMassProj->GetYaxis()->SetRangeUser(10, 150000.0);
    hCorrSideMassRec->Draw("HIST SAME");
    leg3->Draw("SAME");
    
    TCanvas *cCorrSideSideComp = new TCanvas("cCorrSideSideComp", "cCorrSideSideComp", 600, 600);
    cCorrSideSideComp->SetLogy();
    hCorrSideSideProj->Draw("HIST");
    //hCorrSideSideProj->GetYaxis()->SetRangeUser(1, 150000.0);
    hCorrBBProj->Draw("HIST SAME");
    leg4->Draw("SAME");
    
    TCanvas *cRatioMassMass = new TCanvas("cRatioMassMass", "cRatioMassMass", 600, 600);
    hRatioMassMass->Draw("");
    
    TCanvas *cRatioMassSide = new TCanvas("cRatioMassSide", "cRatioMassSide", 600, 600);
    hRatioMassSide->Draw("");
    
    TCanvas *cRatioSideMass = new TCanvas("cRatioSideMass", "cRatioSideMass", 600, 600);
    hRatioSideMass->Draw("");
    
    TCanvas *cRatioSideSide = new TCanvas("cRatioSideSide", "cRatioSideSide", 600, 600);
    hRatioSideSide->Fit("pol0");
    hRatioSideSide->Draw("");
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
