TString sInputNames[2] = {"300m_pThard-1_corr-ntrigg.root", "300m_pThard-1_leading.root"};

void CompareLeadingCorrelations()
{

    gStyle->SetOptStat(0);

    TFile *fIn1 = TFile::Open(sInputNames[0]);
    TFile *fIn2 = TFile::Open(sInputNames[1]);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hCounter = (TH1D*)fIn1->Get("hCounter");
    int nEvent = hCounter->GetBinContent(1);
    std::cout << "Input file contains " << nEvent << " events, proceed to analyse" << std::endl;
    
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
    //hCorrMassMassProj->Scale(1./0.56);
    hCorrMassMassProj->SetTitle("f_{mass,mass} comparison");
            
    TH1D *hCorrMassSideProj = hCorrMassSide->ProjectionX();
    hCorrMassSideProj->SetLineColor(kBlack);
    //hCorrMassSideProj->Scale(1./0.61);
    hCorrMassSideProj->SetTitle("f_{mass,side} comparison");
            
    TH1D *hCorrSideMassProj = hCorrSideMass->ProjectionX();
    hCorrSideMassProj->SetLineColor(kBlack);
    //hCorrSideMassProj->Scale(1./1.8);
    hCorrSideMassProj->SetTitle("f_{side,mass} comparison");
            
    TH1D *hCorrSideSideProj = hCorrSideSide->ProjectionX();
    hCorrSideSideProj->SetLineColor(kBlack);
    //hCorrSideSideProj->Scale(1./4.9);
    hCorrSideSideProj->SetTitle("f_{side,side} comparison");

    // real correlations
    TH1D *hCorrSSProj = hCorrSS->ProjectionX();
    TH1D *hCorrSBProj = hCorrSB->ProjectionX();
    TH1D *hCorrBSProj = hCorrBS->ProjectionX();
    TH1D *hCorrBBProj = hCorrBB->ProjectionX();
    hCorrBBProj->SetLineColor(kRed);
    hCorrBBProj->Scale(4.9);

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
    hCorrMassMassProj->GetYaxis()->SetRangeUser(50, 150000.0);
    leg1->Draw("SAME");
    
    TCanvas *cCorrMassSideComp = new TCanvas("cCorrMassSideComp", "cCorrMassSideComp", 600, 600);
    cCorrMassSideComp->SetLogy();
    hCorrMassSideProj->Draw("HIST");
    hCorrMassSideRec->Draw("HIST SAME");
    hCorrMassSideProj->GetYaxis()->SetRangeUser(10, 150000.0);
    leg2->Draw("SAME");
    
    TCanvas *cCorrSideMassComp = new TCanvas("cCorrSideMassComp", "cCorrSideMassComp", 600, 600);
    cCorrSideMassComp->SetLogy();
    hCorrSideMassProj->Draw("HIST");
    hCorrSideMassProj->GetYaxis()->SetRangeUser(10, 150000.0);
    hCorrSideMassRec->Draw("HIST SAME");
    leg3->Draw("SAME");
    
    TCanvas *cCorrSideSideComp = new TCanvas("cCorrSideSideComp", "cCorrSideSideComp", 600, 600);
    cCorrSideSideComp->SetLogy();
    hCorrSideSideProj->Draw("HIST");
    hCorrSideSideProj->GetYaxis()->SetRangeUser(1, 150000.0);
    hCorrBBProj->Draw("HIST SAME");
    leg4->Draw("SAME");
    
    TCanvas *cRatioMassMass = new TCanvas("cRatioMassMass", "cRatioMassMass", 600, 600);
    hRatioMassMass->Draw("");
    
    TCanvas *cRatioMassSide = new TCanvas("cRatioMassSide", "cRatioMassSide", 600, 600);
    hRatioMassSide->Draw("");
    
    TCanvas *cRatioSideMass = new TCanvas("cRatioSideMass", "cRatioSideMass", 600, 600);
    hRatioSideMass->Draw("");
    
    TCanvas *cRatioSideSide = new TCanvas("cRatioSideSide", "cRatioSideSide", 600, 600);
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
