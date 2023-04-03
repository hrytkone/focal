
const int nset = 3;
const TString filename[nset] = {
    "/home/heimarry/Simulations/focal/analysis_output_2/2023-01-26_pp-focal_no-photon-eff.root",
    "/home/heimarry/Simulations/focal/analysis_output_2/2023-01-27_pp-focal_asym-08.root",
    "/home/heimarry/Simulations/focal/analysis_output_2/2023-01-31_pp-focal_asym-05.root"
};

const int nTriggBins = 2;
const int nAssocBins = 2;

const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

const TString label[nset] = {"no asymmetry cut", Form("#alpha_{thres} = %0.1f", 0.8), Form("#alpha_{thres} = %0.1f", 0.5)};
const EColor mColor[nset] = {kBlack, kBlue, kRed};

TFile *fin[nset];
TH1D *hCounter[nset];
TH1D *hMassAssoc[nset][nTriggBins][nAssocBins];
TCanvas *c1;
TLegend *leg1;

void SetStyle(Bool_t graypalette);
void LoadData();
void ConfigData();
void PlotData();

void PlotMassWithDifferentAlpha()
{
    gStyle->SetOptStat(0);
    SetStyle(0);
    LoadData();
    ConfigData();

    leg1 = new TLegend(0.45, 0.65, 0.75, 0.9);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.038);
    leg1->SetHeader(Form("#splitline{PYTHIA8 simulation}{p_{T} = [%4.1f,%4.1f][%4.1f,%4.1f] GeV/c}",triggPt[0],triggPt[1],assocPt[0],assocPt[1]));
    leg1->AddEntry(hMassAssoc[0][0][0], " ", "");
    for (int i=0; i<nset; i++) leg1->AddEntry(hMassAssoc[i][0][0], label[i], "lep");

    PlotData();
    c1->SaveAs("asymmetry_mass.pdf");
}

//************************************************************************************************
//********************************  FUNCTIONS  ***************************************************
//************************************************************************************************

void LoadData()
{
    for (int i = 0; i < nset; i++) {
        fin[i] = TFile::Open(filename[i]);
        hCounter[i] = (TH1D*)fin[i]->Get("hCounter");
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;
                hMassAssoc[i][itrigg][iassoc] = (TH1D*)fin[i]->Get(Form("Masses/hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            }
        }
    }
}

void ConfigData()
{
    for (int i = 0; i < nset; i++) {
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;
                hMassAssoc[i][itrigg][iassoc]->SetLineColor(mColor[i]);
                hMassAssoc[i][itrigg][iassoc]->SetLineWidth(2);
                hMassAssoc[i][itrigg][iassoc]->Scale(1./hCounter[i]->GetBinContent(1), "width");
                hMassAssoc[i][itrigg][iassoc]->GetXaxis()->SetRangeUser(0.1, 395.);
                hMassAssoc[i][itrigg][iassoc]->GetXaxis()->SetTitleSize(0.045);
                hMassAssoc[i][itrigg][iassoc]->GetXaxis()->SetLabelSize(0.045);
                hMassAssoc[i][itrigg][iassoc]->GetYaxis()->SetTitleSize(0.045);
                hMassAssoc[i][itrigg][iassoc]->GetYaxis()->SetLabelSize(0.045);
                hMassAssoc[i][itrigg][iassoc]->SetTitle(";m_{#gamma#gamma}; 1/N_{ev} dN/dm_{#gamma#gamma}");
                hMassAssoc[i][itrigg][iassoc]->SetMarkerColor(mColor[i]);
                //hMassAssoc[i][itrigg][iassoc]->SetFillColor(mColor[i]);
            }
        }
    }
}

void PlotData()
{
    c1 = new TCanvas("c1", "c1", 600, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.1);
    for (int i = 0; i < nset; i++) {
        if (i==0)
           hMassAssoc[i][0][0]->Draw("E HIST");
        else
           hMassAssoc[i][0][0]->Draw("E HIST SAME");
    }
    leg1->Draw("SAME");
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
