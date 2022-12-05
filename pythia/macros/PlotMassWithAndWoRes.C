const int nset = 2;
const int npt = 2;
const TString filenames[nset] = {"../mass-no-res.root", "../mass-res.root"};
const double triggPt[npt] = {1.0, 10000.0};

double lmargin[nset] = {0.10, 0.12};
double rmargin[nset] = {0.02, 0.02};
double tmargin[nset] = {0.05, 0.05};
double bmargin[nset] = {0.10, 0.10};

TFile *fin[nset];
TH1D *hMass[nset];

void LoadData();
void SetStyle(Bool_t graypalette);

void PlotMassWithAndWoRes()
{
    LoadData();
    SetStyle(0);

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->Divide(2,1);

    for (int i=0; i<nset; i++) {
        c1->cd(i+1);
        gPad->SetLeftMargin(lmargin[i]);
        gPad->SetRightMargin(rmargin[i]);
        gPad->SetTopMargin(tmargin[i]);
        gPad->SetBottomMargin(bmargin[i]);
        hMass[i]->Draw("HIST E");
    }
    c1->SaveAs("mass-example.pdf");
}

//******************************************************************************
//******************************************************************************

void LoadData()
{
    for (int i=0; i<nset; i++) {
        fin[i] = TFile::Open(filenames[i].Data(), "READ");
        cout << Form("hPi0MassTrigg[%4.1f,%4.1f]",triggPt[0],triggPt[1]) << endl;
        hMass[i] = (TH1D*)fin[i]->Get(Form("Masses/hPi0MassTrigg[%4.1f,%4.1f]",triggPt[0],triggPt[1]));
        hMass[i]->Rebin();
        hMass[i]->GetXaxis()->SetRangeUser(0., 390);
        hMass[i]->GetXaxis()->SetLabelSize(0.042);
        hMass[i]->GetYaxis()->SetLabelSize(0.042);
        hMass[i]->GetXaxis()->SetTitleSize(0.048);
        hMass[i]->GetYaxis()->SetTitleSize(0.048);
        hMass[i]->GetYaxis()->SetMaxDigits(3);
        hMass[i]->SetLineColor(kBlack);
        hMass[i]->SetLineWidth(2);
        hMass[i]->SetFillColor(kGray);
        hMass[i]->SetTitle(";m_{#gamma#gamma};counts");
    }
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

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
    //gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    //gStyle->SetFuncColor(kGreen);
    //gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.055,"xyz");
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
