void PlotComponents(TString sInputName = "output.root")
{

    //gStyle->SetOptStat(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hCounter = (TH1D*)fIn->Get("hCounter");
    int nEvent = hCounter->GetBinContent(1);
    std::cout << "Input file contains " << nEvent << " events, proceed to analyse" << std::endl;
    //int nRealTrigg = hCounter->GetBinContent(2);
    //int nMeasTrigg = hCounter->GetBinContent(3);
    //int nRealRecTrigg = hCounter->GetBinContent(4);
    //int nFakeRecTrigg = hCounter->GetBinContent(5);
    //std::cout << "Real triggers : " << nRealTrigg << std::endl;
    //std::cout << "Measured triggers : " << nMeasTrigg << std::endl;
    //std::cout << "Real (rec) triggers : " << nRealRecTrigg << std::endl;
    //std::cout << "Real (fake) triggers : " << nFakeRecTrigg << std::endl;

    TH2D* hCorrMeasured = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",2.,4.,1.5,2.));
    //hCorrMeasured->Rebin2D(4);
    hCorrMeasured->SetLineColor(kBlack);
    hCorrMeasured->GetYaxis()->SetMaxDigits(3);

    TH2D* hCorrSS = (TH2D*)fIn->Get(Form("TrueComponents/hCorrSignalSignal[%4.1f,%4.1f][%4.1f,%4.1f]",2.,4.,1.5,2.));
    //hCorrSS->Rebin2D(4);
    hCorrSS->GetYaxis()->SetMaxDigits(3);

    TH2D* hCorrSB = (TH2D*)fIn->Get(Form("TrueComponents/hCorrSignalBg[%4.1f,%4.1f][%4.1f,%4.1f]",2.,4.,1.5,2.));
    //hCorrSB->Rebin2D(4);
    hCorrSB->GetYaxis()->SetMaxDigits(3);

    TH2D* hCorrBS = (TH2D*)fIn->Get(Form("TrueComponents/hCorrBgSignal[%4.1f,%4.1f][%4.1f,%4.1f]",2.,4.,1.5,2.));
    //hCorrBS->Rebin2D(4);
    hCorrBS->GetYaxis()->SetMaxDigits(3);

    TH2D* hCorrBB = (TH2D*)fIn->Get(Form("TrueComponents/hCorrBgBg[%4.1f,%4.1f][%4.1f,%4.1f]",2.,4.,1.5,2.));
    //hCorrBB->Rebin2D(4);
    hCorrBB->GetYaxis()->SetMaxDigits(3);

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.58, 0.65, 0.78, 0.85);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.05);

    // ----------------
    // |   Analysis   |
    // ----------------

    TCanvas *cCorr = new TCanvas("cCorr", "cCorr", 600, 600);

    TH1D* hCorrMeasuredProj = hCorrMeasured->ProjectionX();
    //hCorrMeasuredProj->Scale(1./nMeasTrigg);
    hCorrMeasuredProj->SetTitle("; #Delta#phi; Counts");
    hCorrMeasuredProj->GetYaxis()->SetTitleOffset(1.);
    //hCorrMeasuredProj->GetYaxis()->SetRangeUser(0.00001, 0.6);
    hCorrMeasuredProj->SetLineColor(kGray);
    hCorrMeasuredProj->SetLineWidth(2);

    TH1D* hCorrSSProj = hCorrSS->ProjectionX();
    hCorrSSProj->SetMarkerStyle(7);
    hCorrSSProj->SetMarkerColor(kBlack);
    hCorrSSProj->SetLineColor(kBlack);
    hCorrSSProj->SetLineWidth(2);
    hCorrSSProj->SetLineStyle(1);
    //hCorrSSProj->Scale(1./nRealRecTrigg);
    hCorrSSProj->SetTitle("; #Delta#phi; Counts");

    TH1D* hCorrSBProj = hCorrSB->ProjectionX();
    hCorrSBProj->SetMarkerStyle(7);
    //hCorrSBProj->Scale(1./nRealRecTrigg);
    hCorrSBProj->SetMarkerColor(kBlue);
    hCorrSBProj->SetLineColor(kBlue);
    hCorrSBProj->SetLineWidth(2);
    hCorrSBProj->SetLineStyle(2);

    TH1D* hCorrBSProj = hCorrBS->ProjectionX();
    hCorrBSProj->SetMarkerStyle(7);
    //hCorrBSProj->Scale(1./nFakeRecTrigg);
    hCorrBSProj->SetMarkerColor(kRed);
    hCorrBSProj->SetLineColor(kRed);
    hCorrBSProj->SetTitle("; #Delta#phi; Counts");
    hCorrBSProj->SetLineWidth(2);
    hCorrBSProj->SetLineStyle(2);

    TH1D* hCorrBBProj = hCorrBB->ProjectionX();
    hCorrBBProj->SetMarkerStyle(7);
    //hCorrBBProj->Scale(1./nFakeRecTrigg);
    hCorrBBProj->SetMarkerColor(8);
    hCorrBBProj->SetLineColor(8);
    hCorrBBProj->SetLineWidth(2);
    hCorrBBProj->SetLineStyle(2);

    leg1->AddEntry(hCorrMeasuredProj, "f_{mass,mass}", "l");
    leg1->AddEntry(hCorrSSProj, "f_{SS}", "l");
    leg1->AddEntry(hCorrSBProj, "f_{SB}", "l");
    leg1->AddEntry(hCorrBSProj, "f_{BS}", "l");
    leg1->AddEntry(hCorrBBProj, "f_{BB}", "l");

    hCorrMeasuredProj->Draw("HIST");
    hCorrSSProj->Draw("HIST SAME");
    hCorrBSProj->Draw("HIST SAME");
    hCorrSBProj->Draw("HIST SAME");
    hCorrBBProj->Draw("HIST SAME");
    leg1->Draw("SAME");

    TCanvas *cBBperSB = new TCanvas("cBBperSB", "cBBperSB", 600, 600);
    TH1D *hBBperSB = (TH1D*)hCorrBBProj->Clone("hBBperSB");
    hBBperSB->SetLineStyle(1);
    hBBperSB->SetLineColor(kBlack);
    hBBperSB->SetMarkerColor(kRed);
    hBBperSB->SetLineWidth(1);
    hBBperSB->Divide(hCorrSBProj);
    hBBperSB->Draw("P");

    TCanvas *cSSperMeas = new TCanvas("cMeasperSS", "cSSperMeas", 600, 600);
    TH1D *hSSperMeas = (TH1D*)hCorrSSProj->Clone("hSSperMeas");
    hSSperMeas->SetLineStyle(1);
    hSSperMeas->SetLineColor(kBlack);
    hSSperMeas->SetMarkerColor(kRed);
    hSSperMeas->SetLineWidth(1);
    hSSperMeas->Divide(hCorrMeasuredProj);
    hSSperMeas->Draw("P");
}


// -----------------
// |   Functions   |
// -----------------
