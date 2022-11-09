const int nset = 2;
const TString input[nset] = {
    "gamma-energy_12900.root",
    "gamma-energy_11220.root"
};
const TString title[nset] = {
    "Calib.Const1=12900",
    "Calib.Const1=11220"
};

TFile *fin[nset];
TH2D *hEnergy[nset];
TCanvas *cEnergy;

void LoadData();
void ConfigPlots();
void Plot();

void PlotGammaEnergyLinearity()
{
    gStyle->SetOptStat(0);
    LoadData();
    ConfigPlots();
    Plot();
}

//______________________________________________________________________________
void LoadData()
{
    for (int iset=0; iset<nset; iset++) {
        fin[iset] = TFile::Open(input[iset].Data());
        hEnergy[iset] = (TH2D*)fin[iset]->Get("hEnergy");
    }
}

void ConfigPlots()
{
    for (int iset=0; iset<nset; iset++) {
        hEnergy[iset]->SetTitle(Form("%s;E_{clust}(MeV);E_{MC}(MeV)", title[iset].Data()));
    }
}

void Plot()
{
    cEnergy = new TCanvas("c1", "", 1200, 600);
    cEnergy->Divide(2,1,0.0001,0.0001);
    TLine *l = new TLine(0, 0, 1600, 1600);
    l->SetLineColor(kRed);
    l->SetLineWidth(2);
    for (int iset=0; iset<nset; iset++) {
        cEnergy->cd(iset+1);
        //gPad->SetLogy();
        //gPad->SetLogx();
        gPad->SetLogz();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        hEnergy[iset]->Draw("COL");
        l->Draw("SAME");
    }
}
