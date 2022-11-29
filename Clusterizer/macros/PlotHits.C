const int nlayers = 20;

TFile *fin;
TH2D* hHitsEnergy[nlayers];
TCanvas *c1[nlayers];

void LoadData(TString input);
void Plot();

void PlotHits(TString input)
{
    gStyle->SetCanvasPreferGL(true);
    gStyle->SetOptStat(0);

    LoadData(input);
    Plot();
}

//******************************************************************************
//******************************************************************************

void LoadData(TString input)
{
    fin = TFile::Open(input.Data(), "READ");
    for (int i=0; i<nlayers; i++)
        hHitsEnergy[i] = (TH2D*)fin->Get(Form("hHitEnergyMap_L%d", i));
}

void Plot()
{
    for (int i=0; i<nlayers; i++) {
        c1[i] = new TCanvas(Form("c1_%d", i), "c1", 600, 600);
        TLatex lpix;
        if (i==4 || i==9) {
            lpix.SetTextAlign(23);
            lpix.SetTextSize(0.08);
            lpix.DrawLatexNDC(0.5, 0.5, "PIXEL");
        }

        hHitsEnergy[i]->Draw("COL");
        c1[i]->SaveAs(Form("hits_%d.png", i));
    }
}
