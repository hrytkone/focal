const int nset = 4;
const int nseg = 6;
const int nwidth = 2;

const TString filename[nset] = {
    "shower-shape_pi0.root",
    "shower-shape_piplus.root",
    "shower-shape_gamma.root",
    "shower-shape_electron.root"
};
const TString legHeader[nseg][nwidth] = {
    {"#splitline{FoCal-E shower shape 1}{Segment 0, LG}", "#splitline{FoCal-E shower shape 2}{Segment 0, LG}"},
    {"#splitline{FoCal-E shower shape 1}{Segment 1, HG}", "#splitline{FoCal-E shower shape 2}{Segment 1, HG}"},
    {"#splitline{FoCal-E shower shape 1}{Segment 2, LG}", "#splitline{FoCal-E shower shape 2}{Segment 2, LG}"},
    {"#splitline{FoCal-E shower shape 1}{Segment 3, HG}", "#splitline{FoCal-E shower shape 2}{Segment 3, HG}"},
    {"#splitline{FoCal-E shower shape 1}{Segment 4, LG}", "#splitline{FoCal-E shower shape 2}{Segment 4, LG}"},
    {"#splitline{FoCal-E shower shape 1}{Segment 5, LG}", "#splitline{FoCal-E shower shape 2}{Segment 5, LG}"}
};
const TString legHeader2[nseg] = {
    "#splitline{FoCal-E shower shape}{Segment 0, LG}",
    "#splitline{FoCal-E shower shape}{Segment 1, HG}",
    "#splitline{FoCal-E shower shape}{Segment 2, LG}",
    "#splitline{FoCal-E shower shape}{Segment 3, HG}",
    "#splitline{FoCal-E shower shape}{Segment 4, LG}",
    "#splitline{FoCal-E shower shape}{Segment 5, LG}"
};
const TString legEntry1[nset] = {
    "#pi^{0}",
    "#pi^{+}",
    "#gamma",
    "e^{-}"
};

const int mColor[nset] = {910, 860, 1, 825};

const double xi[nseg] = {0.42, 0.42, 0.42, 0.42, 0.42, 0.42};
const double xf[nseg] = {0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
const double yi[nseg] = {0.55, 0.55, 0.55, 0.55, 0.55, 0.55};
const double yf[nseg] = {0.87, 0.87, 0.87, 0.87, 0.87, 0.87};

const double xi2[nseg] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
const double xf2[nseg] = {0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
const double yi2[nseg] = {0.55, 0.55, 0.55, 0.55, 0.55, 0.55};
const double yf2[nseg] = {0.87, 0.87, 0.87, 0.87, 0.87, 0.87};

//------------------------------------------------------------------------------

TFile *fin[nset];
TH1D *hShowerShape[nset][nseg][nwidth];
TH1D *hShowerShapeComparison[nset][nseg];
TLegend *leg[nseg][nwidth];
TLegend *leg2[nseg];
TCanvas *c1[nwidth];
TCanvas *c2;

void LoadData();
void CreateLegends();

void PlotShowerShape()
{
    gStyle->SetOptStat(0);

    LoadData();
    CreateLegends();
    for (int k=0; k<nwidth; k++) {
        c1[k] = new TCanvas(Form("c_%d", k), "", 900, 600);
        c1[k]->Divide(3,2);
        for (int j=0; j<nseg; j++) {
            c1[k]->cd(j+1);
            gPad->SetLeftMargin(0.1);
            gPad->SetBottomMargin(0.1);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            for (int i=0; i<nset; i++) {
                hShowerShape[i][j][k]->SetTitle(";w;dN/dw");
                double max = hShowerShape[i][j][k]->GetBinContent(hShowerShape[i][j][k]->GetMaximumBin());
                hShowerShape[i][j][k]->GetYaxis()->SetRangeUser(0., max+0.2*max);
                if (i==0)
                    hShowerShape[i][j][k]->Draw("PE");
                else
                    hShowerShape[i][j][k]->Draw("PE SAME");
                leg[j][k]->Draw("SAME");
            }
        }
    }

    c2 = new TCanvas("c2", "", 900, 600);
    c2->Divide(3,2);
    for (int j=0; j<nseg; j++) {
        c2->cd(j+1);
        gPad->SetLeftMargin(0.1);
        gPad->SetBottomMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        for (int i=0; i<nset; i++) {
            hShowerShapeComparison[i][j]->SetTitle(";w_{2}/w_{1};dN/d(w_{2}/w_{1})");
            double max = hShowerShapeComparison[i][j]->GetBinContent(hShowerShapeComparison[i][j]->GetMaximumBin());
            hShowerShapeComparison[i][j]->GetYaxis()->SetRangeUser(0., max+0.2*max);
            if (i==0)
                hShowerShapeComparison[i][j]->Draw("PE");
            else
                hShowerShapeComparison[i][j]->Draw("PE SAME");
            leg2[j]->Draw("SAME");
        }
    }
}

void LoadData()
{
    for (int i=0; i<nset; i++) {
        fin[i] = TFile::Open(filename[i].Data());
        for (int j=0; j<nseg; j++) {
            hShowerShapeComparison[i][j] = (TH1D*)fin[i]->Get(Form("hWidth2PerWidth1_%d", j));
            hShowerShapeComparison[i][j]->Rebin(2);
            hShowerShapeComparison[i][j]->SetMarkerStyle(20);
            hShowerShapeComparison[i][j]->SetMarkerSize(0.5);
            hShowerShapeComparison[i][j]->SetMarkerColor(mColor[i]);
            hShowerShapeComparison[i][j]->SetLineColor(mColor[i]);
            hShowerShapeComparison[i][j]->GetYaxis()->SetMaxDigits(3);
            hShowerShapeComparison[i][j]->Scale(1./hShowerShapeComparison[i][j]->GetEntries(), "width");
            for (int k=0; k<nwidth; k++) {
                hShowerShape[i][j][k] = (TH1D*)fin[i]->Get(Form("hWidth%d_%d", k+1, j));
                hShowerShape[i][j][k]->Rebin(2);
                hShowerShape[i][j][k]->SetMarkerStyle(20);
                hShowerShape[i][j][k]->SetMarkerSize(0.5);
                hShowerShape[i][j][k]->SetMarkerColor(mColor[i]);
                hShowerShape[i][j][k]->SetLineColor(mColor[i]);
                hShowerShape[i][j][k]->GetYaxis()->SetMaxDigits(3);
                hShowerShape[i][j][k]->Scale(1./hShowerShape[i][j][k]->GetEntries(), "width");
            }
        }
    }
}

void CreateLegends()
{
    for (int iseg=0; iseg<nseg; iseg++) {
        for (int k=0; k<nwidth; k++) {
            leg[iseg][k] = new TLegend(xi[iseg], yi[iseg], xf[iseg], yf[iseg]);
            leg[iseg][k]->SetFillStyle(0); leg[iseg][k]->SetBorderSize(0); leg[iseg][k]->SetTextSize(gStyle->GetTextSize()*1.);
            leg[iseg][k]->SetHeader(legHeader[iseg][k]);
            leg[iseg][k]->AddEntry(hShowerShape[0][iseg][k], "", "");
            for (int iset=0; iset<nset; iset++)
                leg[iseg][k]->AddEntry(hShowerShape[iset][iseg][k], legEntry1[iset], "ep");
            TLegendEntry *header = (TLegendEntry*)leg[iseg][k]->GetListOfPrimitives()->First();
            header->SetTextAlign(11);
        }
        leg2[iseg] = new TLegend(xi2[iseg], yi2[iseg], xf2[iseg], yf2[iseg]);
        leg2[iseg]->SetFillStyle(0); leg2[iseg]->SetBorderSize(0); leg2[iseg]->SetTextSize(gStyle->GetTextSize()*1.);
        leg2[iseg]->SetHeader(legHeader2[iseg]);
        leg2[iseg]->AddEntry(hShowerShape[0][iseg][0], "", "");
        for (int iset=0; iset<nset; iset++)
            leg2[iseg]->AddEntry(hShowerShapeComparison[iset][iseg], legEntry1[iset], "ep");
        TLegendEntry *header = (TLegendEntry*)leg2[iseg]->GetListOfPrimitives()->First();
        header->SetTextAlign(11);
    }
}
