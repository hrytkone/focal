const int netabins = 13;

void SetStyle(Bool_t graypalette);

void PlotEffEtaBins(TString sInputName = "input.root")
{
    gStyle->SetOptStat(0);
    SetStyle(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hPionPtDetected[netabins];
    TH1D *hPionPtFor[netabins];
    TH1D *hRatio[netabins];
    for (int i=0; i<netabins; i++) {
        hPionPtDetected[i] = (TH1D*)fIn->Get(Form("hPionPtDetected_%d", i));
        hPionPtFor[i]      = (TH1D*)fIn->Get(Form("hPionPtFor_%d", i));
        hRatio[i] = (TH1D*)hPionPtDetected[i]->Clone(Form("hRatio_%d", i));
        hRatio[i]->Divide(hPionPtFor[i]);
    }

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.28, 0.25, 0.48, 0.45);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.035);
    leg1->SetHeader("#eta bins (bin width = 0.2)");
    leg1->AddEntry(hRatio[0], "[3.2, 3.4]", "pe");
    leg1->AddEntry(hRatio[0], "#upoint#upoint#upoint", "");
    leg1->AddEntry(hRatio[netabins-1], "[5.6, 5.8]", "pe");

    // ----------------
    // |   Analysis   |
    // ----------------
    hRatio[0]->SetTitle(" ; p_{T} (GeV); Pair efficiency");
    hRatio[0]->GetYaxis()->SetTitleOffset(1.);
    hRatio[0]->GetYaxis()->SetRangeUser(0., 1.1);
    //hRatio[0]->GetXaxis()->SetRangeUser(0., 12.);

    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 600, 600);
    cRatio->SetLogx();
    for (int i=0; i<netabins; i++) {
        hRatio[i]->SetMarkerSize(0.5);
        hRatio[i]->SetMarkerStyle(20);
        if (i==0)
            hRatio[i]->Draw("PLC PMC");
        else
            hRatio[i]->Draw("PLC PMC SAME");
    }
    leg1->Draw("SAME");
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    gStyle->SetPalette(kCool);
    //gStyle->Reset("Plain");
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineScalePS(1);
    //if(graypalette) gStyle->SetPalette(8,0);
    //else gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    //gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    //gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.035,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    //gStyle->SetTitleSize(0.035,"xyz");
    //gStyle->SetTitleOffset(1.25,"y");
    //gStyle->SetTitleOffset(1.2,"x");
    //gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    //  gStyle->SetFillColor(kWhite);
    gStyle->SetLegendFont(42);
}
