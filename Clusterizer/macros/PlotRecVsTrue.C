TFile *fin;
TH2D *hPt, *hEta, *hE;

void redrawBorder();
void SetStyle(Bool_t graypalette);

void PlotRecVsTrue(TString input)
{
    SetStyle(0);

    fin = TFile::Open(input.Data(), "READ");
    hPt = (TH2D*)fin->Get("hRecPtVsTruePt");
    hEta = (TH2D*)fin->Get("hRecEtaVsTrueEta");
    hE = (TH2D*)fin->Get("hRecEVsTrueE");

    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 600);
    c1->Divide(3,1);

    c1->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    hPt->Draw("COLZ");
    TLine *l1 = new TLine(2., 2., 20., 20.);
    l1->SetLineColor(kRed);
    l1->SetLineWidth(2);
    l1->Draw("SAME");
    
    c1->cd(2);
    hEta->Draw("COLZ");
    TLine *l2 = new TLine(3., 3., 6., 6.);
    l2->SetLineColor(kRed);
    l2->SetLineWidth(2);
    l2->Draw("SAME");
    
    c1->cd(3);
    gPad->SetLogx();
    gPad->SetLogy();
    hE->Draw("COLZ");
    hE->GetXaxis()->SetRangeUser(20, 2000.);
    hE->GetYaxis()->SetRangeUser(20, 2000.);
    TLine *l3 = new TLine(10., 10., 2000., 2000.);
    l3->SetLineColor(kRed);
    l3->SetLineWidth(2);
    l3->Draw("SAME");
}

void redrawBorder()
{
    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.SetLineWidth(2);
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetLineScalePS(1);
    //gStyle->SetPalette(kCool);
    if(graypalette) gStyle->SetPalette(8,0);
    //else gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    //gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    //gStyle->SetPadTickX(0);
    //gStyle->SetPadTickY(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(1);
    gStyle->SetFuncColor(kRed);
    gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.042,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.042,"xyz");
    gStyle->SetTitleOffset(1.2,"y");
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleFillColor(kWhite);
    //gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    //gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    //  gStyle->SetFillColor(kWhite);
    gStyle->SetLegendFont(42);
}
