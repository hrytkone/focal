const TString inputname = "asymmetry_output.root";

const int nset = 3;
const TString houtname[nset] = {
    "hAsymTrue",
    "hAsymRecPythia",
    "hAsymRecGeant"
};
const TString legentries[nset] = {
    Form("MC truth (#times65)"),
    Form("PYTHIA8 photons (#times20)"),
    Form("PYTHIA6+GEANT3 clusters")
};
TH1D *hist[nset];

const EColor lColor[nset] = {kBlack, kBlack, kRed};
const int lStyle[nset] = {2, 1, 1};
double scale[nset] = {65., 20., 1.0};

void redrawBorder();
void SetStyle(Bool_t graypalette);

void PlotAsymmetry()
{
    SetStyle(0);

    TFile *fin = TFile::Open(inputname.Data());
    for (int i = 0; i < nset; i++) {
        hist[i] = (TH1D*)fin->Get(houtname[i].Data());
        hist[i]->SetTitle("; Asymmetry #alpha; 1/N_{ev} dN/d#alpha");
        //hist[i]->GetXaxis()->SetRangeUser(0.1, 1.);
        hist[i]->SetLineColor(lColor[i]);
        hist[i]->SetMarkerColor(lColor[i]);
        hist[i]->SetLineStyle(lStyle[i]);
        hist[i]->SetLineWidth(2);
        hist[i]->GetXaxis()->SetTitleSize(0.045);
        hist[i]->GetYaxis()->SetTitleSize(0.045);
        hist[i]->GetXaxis()->SetLabelSize(0.045);
        hist[i]->GetYaxis()->SetLabelSize(0.045);
        hist[i]->GetYaxis()->SetTitleOffset(1.4);
        hist[i]->Scale(scale[i], "width");
    }

    TLegend *leg = new TLegend(0.2, 0.2, 0.35, 0.35);
    leg->SetBorderSize(0); leg->SetTextSize(0.042);
    for (int i = 0; i < nset; i++) leg->AddEntry(hist[i], legentries[i], "le");

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.05);
    hist[nset-1]->Draw("HIST E");
    for (int i = 0; i < nset; i++) hist[i]->Draw("HIST E SAME");
    leg->Draw("SAME");

    TLine *line = new TLine(0.8, 0., 0.8, 250.);
    line->SetLineStyle(9);
    line->Draw("SAME");

    redrawBorder();
    c1->SaveAs("asymmetry.pdf");
}

void redrawBorder()
{
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
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
    gStyle->SetTitleOffset(0.75,"y");
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