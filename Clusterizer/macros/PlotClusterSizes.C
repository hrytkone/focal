const int npt = 6;
const int necut = 6;
//const TString filename = "numberOfClusters.root";
const TString filename = "numberOfClusters_lclust-ecut-3.root";
//const TString outputdir = "20220510_gamma-cluster-results";
const TString outputdir = "20220510_gamma-cluster-results_leading-ecut-3";

const TString legEntry[npt] = {
    "2 GeV < p_{T} < 3 GeV",
    "3 GeV < p_{T} < 4 GeV",
    "4 GeV < p_{T} < 8 GeV",
    "8 GeV < p_{T} < 10 GeV",
    "10 GeV < p_{T} < 15 GeV",
    "15 GeV < p_{T} < 20 GeV"
};

const TString legHeader[necut] = {"#gamma particle gun, no cut", "#gamma particle gun, E > 100 MeV",
                                "#gamma particle gun, E > 350 MeV", "#gamma particle gun, E > 500 MeV",
                                "#gamma particle gun, E > 1 GeV", "#gamma particle gun, E > 3 GeV"};
const double pt[npt+1] = {2., 3., 4., 8., 10., 15., 20.};

//------------------------------------------------------------------------------

TFile *fin;
TH1D *hNumberOfClusters[npt][necut];
TH1 *hCumulative[npt][necut];
TH2D *hClusterPerGammaEnergy;
TH1D *hClusterPerGammaEnergy_py[npt];
TH2D *hLeadingClusterDeltaPhiDeltaEta;
TH1D *hClusterPt;
TH1D *hClusterMoreThanOne;

TLegend *leg[necut];
TCanvas *cClusterDensity[necut];
TCanvas *cCumulative[necut];
TCanvas *cClusterRatio;
TCanvas *cClusterPerGammaEnergy;
TCanvas *cClusterPerGammaEnergyProj;
TCanvas *cLeadingClusterDeltaPhiDeltaEta;
EColor mColor[npt] = {kBlue, kRed, kOrange, kBlack, kMagenta, kCyan};

double xmax[necut] = {25., 20., 20., 20., 20., 15.};

void LoadData();
void CreateLegend();
void SetStyle(Bool_t graypalette);

void PlotClusterSizes()
{
    //gStyle->SetOptStat(0);
    //gStyle->SetPadBottomMargin(0.15);
    //gStyle->SetPadLeftMargin(0.15);

    SetStyle(0);

    LoadData();
    CreateLegend();

    for (int iecut=0; iecut<necut; iecut++) {
        cClusterDensity[iecut] = new TCanvas(Form("cClusterDensity_%d", iecut), "", 600, 600);
        cClusterDensity[iecut]->SetLogy();
        for (int i=0; i<npt; i++) {
            hNumberOfClusters[i][iecut]->SetTitle(";N_{clust};1/N_{clust}dN/dN_{clust}");
            hNumberOfClusters[i][iecut]->GetYaxis()->SetTitleOffset(1.);
            //hNumberOfClusters[i][iecut]->GetYaxis()->SetRangeUser(0., 3e3);
            if (i==0) {
                hNumberOfClusters[i][iecut]->Draw("PE");
            } else {
                hNumberOfClusters[i][iecut]->Draw("SAME");
            }
        }
        leg[iecut]->Draw("SAME");
        cClusterDensity[iecut]->SaveAs(Form("%s/numberofclusters_%d.eps", outputdir.Data(), iecut));

        cCumulative[iecut] = new TCanvas(Form("cCumulative_%d", iecut), "", 600, 600);
        for (int i=0; i<npt; i++) {
            if (i==0) {
                hCumulative[i][iecut]->Draw("HIST");
            } else {
                hCumulative[i][iecut]->Draw("HIST SAME");
            }
        }
        leg[iecut]->Draw("SAME");
        cCumulative[iecut]->SaveAs(Form("%s/cumulative_%d.eps", outputdir.Data(), iecut));
    }

    cClusterRatio = new TCanvas("cClusterRatio", "", 600, 600);
    hClusterMoreThanOne->SetTitle(";p_{T};N_{clust>1}/N_{clust,all}");
    hClusterMoreThanOne->GetYaxis()->SetTitleOffset(1.8);
    hClusterMoreThanOne->Draw("PE");
    cClusterRatio->SaveAs(Form("%s/clusterratio.eps", outputdir.Data()));

    cClusterPerGammaEnergy = new TCanvas("cClusterPerGammaEnergy", "", 600, 600);
    cClusterPerGammaEnergy->SetLogx();
    hClusterPerGammaEnergy->SetTitle(";p_{T};E_{clust}/E_{#gamma}");
    hClusterPerGammaEnergy->GetYaxis()->SetTitleOffset(1.8);
    hClusterPerGammaEnergy->Draw("COLZ");
    cClusterPerGammaEnergy->SaveAs(Form("%s/clusterpergamma.eps", outputdir.Data()));

    cClusterPerGammaEnergyProj = new TCanvas("cClusterPerGammaEnergyProj", "", 600, 600);
    for (int i=0; i<npt; i++) {
        hClusterPerGammaEnergy_py[i]->SetTitle(";E_{clust}/E_{#gamma}; 1/N dN/d(E_{clust}/E_{#gamma})");
        hClusterPerGammaEnergy_py[i]->GetYaxis()->SetTitleOffset(1.8);
        if (i==0)
            hClusterPerGammaEnergy_py[i]->Draw("P");
        else
            hClusterPerGammaEnergy_py[i]->Draw("P SAME");
    }
    leg[necut-1]->Draw("SAME");
    cClusterPerGammaEnergyProj->SaveAs(Form("%s/clusterpergamma_py.eps", outputdir.Data()));

    cLeadingClusterDeltaPhiDeltaEta = new TCanvas("cLeadingClusterDeltaPhiDeltaEta", "", 600, 600);
    hLeadingClusterDeltaPhiDeltaEta->SetTitle(";#Delta#phi;#Delta#eta");
    hLeadingClusterDeltaPhiDeltaEta->GetYaxis()->SetTitleOffset(1.8);
    hLeadingClusterDeltaPhiDeltaEta->Draw("COLZ");
    cLeadingClusterDeltaPhiDeltaEta->SaveAs(Form("%s/deltaphieta.eps", outputdir.Data()));
}

void LoadData()
{
    fin = TFile::Open(filename.Data());
    for (int ipt=0; ipt<npt; ipt++) {
        for (int iecut=0; iecut<necut; iecut++) {
            hNumberOfClusters[ipt][iecut] = (TH1D*)fin->Get(Form("hNumberOfClusters_%d_%d", ipt, iecut));
            hNumberOfClusters[ipt][iecut]->Scale(1./hNumberOfClusters[ipt][iecut]->GetEntries());
            hNumberOfClusters[ipt][iecut]->SetMarkerStyle(kOpenCircle);
            hNumberOfClusters[ipt][iecut]->SetMarkerSize(.5);
            hNumberOfClusters[ipt][iecut]->SetMarkerColor(mColor[ipt]);
            hNumberOfClusters[ipt][iecut]->SetLineColor(mColor[ipt]);
            hNumberOfClusters[ipt][iecut]->GetYaxis()->SetMaxDigits(3);
            hNumberOfClusters[ipt][iecut]->SetLineWidth(1);

            hCumulative[ipt][iecut] = hNumberOfClusters[ipt][iecut]->GetCumulative();
            hCumulative[ipt][iecut]->SetMarkerStyle(kOpenCircle);
            hCumulative[ipt][iecut]->SetMarkerSize(.5);
            hCumulative[ipt][iecut]->SetMarkerColor(mColor[ipt]);
            hCumulative[ipt][iecut]->SetLineColor(mColor[ipt]);
            hCumulative[ipt][iecut]->GetYaxis()->SetMaxDigits(3);
            hCumulative[ipt][iecut]->SetTitle(";N_{clust};CDF");
            hCumulative[ipt][iecut]->SetLineWidth(1);
            hCumulative[ipt][iecut]->GetXaxis()->SetRangeUser(0., xmax[iecut]);
            hCumulative[ipt][iecut]->GetYaxis()->SetTitleOffset(1.);
            hCumulative[ipt][iecut]->GetYaxis()->SetTitleSize(0.042);
            hCumulative[ipt][iecut]->GetXaxis()->SetTitleSize(0.042);
            hCumulative[ipt][iecut]->GetYaxis()->SetLabelSize(0.042);
            hCumulative[ipt][iecut]->GetXaxis()->SetLabelSize(0.042);
        }
    }
    hClusterPt = (TH1D*)fin->Get("hClusterPt");
    hClusterMoreThanOne = (TH1D*)fin->Get("hClusterMoreThanOne");
    hClusterMoreThanOne->Divide(hClusterPt);

    hClusterPerGammaEnergy = (TH2D*)fin->Get("hClusterPerGammaEnergy");
    hClusterPerGammaEnergy->Scale(1., "width");

    for (int i=0; i<npt; i++) {
        int minbin = hClusterPerGammaEnergy->GetXaxis()->FindBin(pt[i]);
        int maxbin = hClusterPerGammaEnergy->GetXaxis()->FindBin(pt[i+1]);
        hClusterPerGammaEnergy_py[i] = hClusterPerGammaEnergy->ProjectionY(Form("_py%d", i), minbin, maxbin);
        hClusterPerGammaEnergy_py[i]->Rebin();
        hClusterPerGammaEnergy_py[i]->Scale(1./hClusterPerGammaEnergy_py[i]->Integral());
        cout << "minbin:" << minbin << " maxbin:" << maxbin << " int:" << hClusterPerGammaEnergy_py[i]->Integral() << endl;
        hClusterPerGammaEnergy_py[i]->GetYaxis()->SetRangeUser(0.00001, 0.15);
        hClusterPerGammaEnergy_py[i]->SetMarkerStyle(kOpenCircle);
        hClusterPerGammaEnergy_py[i]->SetMarkerSize(.5);
        hClusterPerGammaEnergy_py[i]->SetMarkerColor(mColor[i]);
        hClusterPerGammaEnergy_py[i]->SetLineColor(mColor[i]);
        hClusterPerGammaEnergy_py[i]->GetYaxis()->SetMaxDigits(3);
    }

    hLeadingClusterDeltaPhiDeltaEta = (TH2D*)fin->Get("hLeadingClusterDeltaPhiDeltaEta");
    //hLeadingClusterDeltaPhiDeltaEta->Scale(1., "width");
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    //gStyle->Reset("Plain");
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetLineScalePS(1);
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
    gStyle->SetHistLineWidth(2);
    //gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    //gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.045,"xyz");
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

void CreateLegend()
{
    // Legend for singal-to-background figure
    for (int iecut=0; iecut<necut; iecut++) {
        leg[iecut] = new TLegend(0.47, 0.27, 0.82, 0.67);
        leg[iecut]->SetFillStyle(0); leg[iecut]->SetBorderSize(0); leg[iecut]->SetTextSize(0.035);
        //leg[iecut]->SetNColumns(2);
        leg[iecut]->SetHeader(legHeader[iecut].Data(), "C");
        for (int i=0; i<npt; i++)
            leg[iecut]->AddEntry(hNumberOfClusters[i][iecut], legEntry[i], "l");
    }
}
