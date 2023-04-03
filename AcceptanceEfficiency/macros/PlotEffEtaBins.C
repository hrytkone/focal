const int netabins = 6;
const double eta[netabins+1] = {3.2, 3.5, 4.0, 4.5, 5.0, 5.5, 5.8};

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
        hPionPtDetected[i]->Rebin();
        hPionPtFor[i] = (TH1D*)fIn->Get(Form("hPionPtFor_%d", i));
        hPionPtFor[i]->Rebin();
        hRatio[i] = (TH1D*)hPionPtDetected[i]->Clone(Form("hRatio_%d", i));
        hRatio[i]->Divide(hPionPtFor[i]);
    }

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.5, 0.2, 0.9, 0.4);
    leg1->SetNColumns(2);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.035);
    leg1->SetHeader("pseudorapidity bins");
    for (int i=0; i<netabins; i++)
        leg1->AddEntry(hRatio[i], Form("[%0.1f, %0.1f]", eta[i], eta[i+1]), "pe");


    // ----------------
    // |   Analysis   |
    // ----------------
    hRatio[0]->SetTitle(" ; p_{T} (GeV); Pair efficiency");
    hRatio[0]->GetXaxis()->SetTitleOffset(1.2);
    hRatio[0]->GetYaxis()->SetTitleOffset(1.2);
    hRatio[0]->GetYaxis()->SetRangeUser(0., 1.1);
    hRatio[0]->GetXaxis()->SetRangeUser(0.1, 10.);
    hRatio[0]->GetXaxis()->SetTitleSize(0.048);
    hRatio[0]->GetYaxis()->SetTitleSize(0.048);
    hRatio[0]->GetXaxis()->SetLabelSize(0.048);
    hRatio[0]->GetYaxis()->SetLabelSize(0.048);

    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 600, 600);
    cRatio->SetLogx();
    for (int i=0; i<netabins; i++) {
        hRatio[i]->SetMarkerSize(1.);
        hRatio[i]->SetLineWidth(2);
        hRatio[i]->SetMarkerStyle(20);
        if (i==0)
            hRatio[i]->Draw("PLC E1 PMC");
        else
            hRatio[i]->Draw("PLC E1 PMC SAME");
    }
    leg1->Draw("SAME");
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    Int_t palette[6];
    palette[0] = 1;
    palette[1] = 600;
    palette[2] = 616;
    palette[3] = 879;
    palette[4] = 632;
    palette[5] = 807;
    gStyle->SetPalette(6,palette);

    //gStyle->Reset("Plain");
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineScalePS(1);
    //if(graypalette) gStyle->SetPalette(8,0);
    //else gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.05);
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
