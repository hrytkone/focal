const int necut = 6;
double ecut[necut] = {0., 0.1, 0.35, 1., 3., 5.};
int mColor[necut] = {433, 38, 38, 596, 590, 603};

const int nset = 1;
const TString input[nset] = {
    "numberOfClusters.root",
    //"numberOfClusters_eta-35-4.root",
    //"numberOfClusters_eta-4-45.root",
    //"numberOfClusters_eta-45-5.root",
    //"numberOfClusters_eta-5-55.root"
};

TFile *fin;
TH1D *hNumberOfClusters[nset][necut];
TLegend *leg;

void LoadData();
void ConfigHistos();
void CalculatePercentages();
void CreateLegend();
void PlotData();
void redrawBorder();
void SetStyle(Bool_t graypalette);

void PlotNumberOfClusters()
{
    SetStyle(0);
    LoadData();
    CreateLegend();
    CalculatePercentages();    
    ConfigHistos();
    PlotData();
}

//******************************************************************************
//******************************************************************************
void LoadData()
{
    for (int iset=0; iset < nset; iset++) {
        fin = TFile::Open(input[iset].Data(), "READ");
        for (int iecut=0; iecut < necut; iecut++) {
            hNumberOfClusters[iset][iecut] = (TH1D*)fin->Get(Form("hNumberOfClustersAll_%d", iecut));
        }
    }
}

void ConfigHistos()
{
    for (int iset=0; iset<nset; iset++) {
        for (int iecut=0; iecut<necut; iecut++) {
            hNumberOfClusters[iset][iecut]->Scale(1./hNumberOfClusters[iset][iecut]->GetEntries());
            hNumberOfClusters[iset][iecut]->SetFillColor(mColor[iecut]);
            hNumberOfClusters[iset][iecut]->SetMarkerSize(0.);
            hNumberOfClusters[iset][iecut]->SetLineWidth(0.);
            //hNumberOfClusters[iset][iecut]->SetLineColor(mColor[iecut]);
            //hNumberOfClusters[iset][iecut]->SetMarkerColor(mColor[iecut]);
            hNumberOfClusters[iset][iecut]->SetMarkerStyle(kFullCircle);
            //hNumberOfClusters[iset][iecut]->SetLineWidth(2);
            hNumberOfClusters[iset][iecut]->GetYaxis()->SetMaxDigits(3);
            hNumberOfClusters[iset][iecut]->GetXaxis()->SetRangeUser(-0.5, 5.5);
            hNumberOfClusters[iset][iecut]->GetYaxis()->SetRangeUser(0., 1.);
            hNumberOfClusters[iset][iecut]->SetTitle(";N_{cluster};N_{cluster}/N_{cluster,all}");
            hNumberOfClusters[iset][iecut]->GetXaxis()->SetTitleSize(0.06);
            hNumberOfClusters[iset][iecut]->GetYaxis()->SetTitleSize(0.06);
            hNumberOfClusters[iset][iecut]->GetXaxis()->SetLabelSize(0.052);
            hNumberOfClusters[iset][iecut]->GetYaxis()->SetLabelSize(0.052);
            hNumberOfClusters[iset][iecut]->GetXaxis()->SetTitleOffset(0.82);
            hNumberOfClusters[iset][iecut]->GetYaxis()->SetTitleOffset(0.82);
        }
    }
}

void CalculatePercentages()
{
    int nbin = 15;
    for (int iset=0; iset<nset; iset++) {
        for (int iecut=0; iecut<necut; iecut++) {
            cout << "\nE_cut = " << ecut[iecut] << endl;
            int allClusters = hNumberOfClusters[iset][iecut]->GetEntries();
            for (int ibin = 1; ibin <= nbin; ibin++) {
                cout << "\tnclust=" << ibin-1 << "\t" << (double)hNumberOfClusters[iset][iecut]->GetBinContent(ibin)/allClusters << endl;
            }
        }
    }
}

void CreateLegend()
{
    leg = new TLegend(0.5, 0.5, 0.85, 0.85);
    leg->SetLineStyle(0); leg->SetTextSize(0.048);
    leg->SetHeader("#splitline{5.4 < #eta_{true} < 5.8}{2 GeV/c < p_{T,true} < 20 GeV/c}");
    leg->AddEntry(hNumberOfClusters[0][0], " ", "");
    leg->AddEntry(hNumberOfClusters[0][0], "No energy cut", "f");
    leg->AddEntry(hNumberOfClusters[0][3], "Energy > 1 GeV", "f");
    leg->AddEntry(hNumberOfClusters[0][4], "Energy > 3 GeV", "f");
    leg->AddEntry(hNumberOfClusters[0][5], "Energy > 5 GeV", "f");
}

void PlotData()
{
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //c1->SetLogy();
   // c1->SetGrayscale();
    hNumberOfClusters[0][0]->SetBarWidth(0.18);
    hNumberOfClusters[0][0]->SetBarOffset(0.1);

    hNumberOfClusters[0][3]->SetBarWidth(0.18);
    hNumberOfClusters[0][3]->SetBarOffset(0.3);

    hNumberOfClusters[0][4]->SetBarWidth(0.18);
    hNumberOfClusters[0][4]->SetBarOffset(0.5);
    
    hNumberOfClusters[0][5]->SetBarWidth(0.18);
    hNumberOfClusters[0][5]->SetBarOffset(0.7);
    
    hNumberOfClusters[0][0]->Draw("BAR");
    hNumberOfClusters[0][3]->Draw("BAR SAME");
    hNumberOfClusters[0][4]->Draw("BAR SAME");
    hNumberOfClusters[0][5]->Draw("BAR SAME");
    leg->Draw("SAME");
    redrawBorder();

}

void redrawBorder()
{
    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.SetLineWidth(2);
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymin());
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
