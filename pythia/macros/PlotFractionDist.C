const TString input = "/home/heimarry/Simulations/focal/analysis_output_2/2023-03-28_pp-focal_asym-08.root";

const int nTriggBins = 2;
const int nAssocBins = 2;
const double triggPt[nTriggBins+1] = {4.0, 8.0, 16.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};
const EColor mColor[nTriggBins][nAssocBins] = {{kBlack, kRed}, {kBlack, kRed}};
const int mStyle[nTriggBins][nAssocBins] = {{71, 71}, {90, 90}};

TFile *fin;
TH1D *hX1[nTriggBins][nAssocBins], *hX2[nTriggBins][nAssocBins], *hX2copy[nTriggBins][nAssocBins];
TH2D *hX1X2[nTriggBins][nAssocBins];
TLegend *leg, *leg2;

double xi[nTriggBins][nAssocBins] = {{0.18, 0.07}, {0.18,0.07}};
double yi[nTriggBins][nAssocBins] = {{0.75, 0.75}, {0.85,0.85}};

void LoadData();
void SetStyle(Bool_t graypalette);
void redrawBorder();

void PlotFractionDist()
{
    SetStyle(0);
    LoadData();

    leg = new TLegend(0.15, 0.4, 0.5, 0.6);
    leg->SetLineWidth(0); leg->SetTextSize(0.062);
    //leg->SetNColumns(2);
    leg->AddEntry(hX1[0][0], "Fraction x_{1}", "pe");
    leg->AddEntry(hX2[0][0], "Fraction x_{2}", "pe");

    leg2 = new TLegend(0.35, 0.2, 0.75, 0.5);
    leg2->SetLineWidth(0); leg2->SetTextSize(0.042);
    leg2->SetHeader("Fraction x_{2}, [p_{T,t}][p_{T,a}]");
    leg2->AddEntry(hX2copy[0][0], Form("[%4.1f,%4.1f][%4.1f,%4.1f] GeV/c",triggPt[0], triggPt[1], assocPt[0], assocPt[1]), "pe");
    leg2->AddEntry(hX2copy[0][1], Form("[%4.1f,%4.1f][%4.1f,%4.1f] GeV/c",triggPt[0], triggPt[1], assocPt[1], assocPt[2]), "pe");
    leg2->AddEntry(hX2copy[1][0], Form("[%4.1f,%4.1f][%4.1f,%4.1f] GeV/c",triggPt[1], triggPt[2], assocPt[0], assocPt[1]), "pe");
    leg2->AddEntry(hX2copy[1][1], Form("[%4.1f,%4.1f][%4.1f,%4.1f] GeV/c",triggPt[1], triggPt[2], assocPt[1], assocPt[2]), "pe");

    TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
    c1->Divide(nTriggBins, nAssocBins);
    int ipanel = 1;
    for (int it=0; it<nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];
        for (int ia=0; ia<nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;
            c1->cd(ipanel);
            gPad->SetLogx();
            if (ipanel==1 || ipanel==3) {
                gPad->SetRightMargin(0.01);
                gPad->SetLeftMargin(0.12);
                hX1[it][ia]->GetXaxis()->SetTitleSize(0.);
            }
            if (ipanel==2 || ipanel==4) {
                gPad->SetLeftMargin(0.01);
                gPad->SetRightMargin(0.12);
                hX1[it][ia]->GetYaxis()->SetLabelSize(0.);
            }
            if (ipanel==1 || ipanel==2) {
                hX1[it][ia]->GetXaxis()->SetTitleSize(0.);
                hX1[it][ia]->GetXaxis()->SetLabelSize(0.);
                gPad->SetBottomMargin(0.012);
                gPad->SetTopMargin(0.12);
            }
            if (ipanel==3 || ipanel==4) {
                hX1[it][ia]->GetYaxis()->SetTitleSize(0.);
                gPad->SetTopMargin(0.005);
                gPad->SetBottomMargin(0.15);
            }


            hX1[it][ia]->Draw("P");
            hX2[it][ia]->Draw("P SAME");
            if (it==0 && ia==0) leg->Draw("SAME");

            TLatex l;
            l.SetTextSize(0.072);
            l.DrawLatexNDC(xi[it][ia], yi[it][ia], Form("p_{T,t} = [%4.1f,%4.1f] GeV/c",tlow,tupp));
            l.DrawLatexNDC(xi[it][ia], yi[it][ia]-0.1, Form("p_{T,a} = [%4.1f,%4.1f] GeV/c",alow,aupp));
            redrawBorder();
            ipanel++;
        }
    }

    TCanvas *c2 = new TCanvas("c2", "", 1200, 800);
    c2->Divide(nTriggBins, nAssocBins);
    ipanel = 1;
    for (int it=0; it<nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];
        for (int ia=0; ia<nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;
            c2->cd(ipanel);
            gPad->SetLogx();
            gPad->SetLogy();
            if (ipanel==1 || ipanel==3) {
                gPad->SetRightMargin(0.01);
                gPad->SetLeftMargin(0.12);
                hX1X2[it][ia]->GetXaxis()->SetTitleSize(0.);
            }
            if (ipanel==2 || ipanel==4) {
                gPad->SetLeftMargin(0.01);
                gPad->SetRightMargin(0.12);
                hX1X2[it][ia]->GetYaxis()->SetLabelSize(0.);
            }
            if (ipanel==1 || ipanel==2) {
                hX1X2[it][ia]->GetXaxis()->SetTitleSize(0.);
                hX1X2[it][ia]->GetXaxis()->SetLabelSize(0.);
                gPad->SetBottomMargin(0.02);
                gPad->SetTopMargin(0.12);
            }
            if (ipanel==3 || ipanel==4) {
                hX1X2[it][ia]->GetYaxis()->SetTitleSize(0.);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.15);
            }

            hX1X2[it][ia]->Draw("COL");

            TLatex l;
            l.SetTextSize(0.082);
            l.DrawLatexNDC(xi[it][ia], yi[it][ia]-0.52, Form("p_{T,t} = [%4.1f,%4.1f] GeV/c",tlow,tupp));
            l.DrawLatexNDC(xi[it][ia], yi[it][ia]-0.62, Form("p_{T,a} = [%4.1f,%4.1f] GeV/c",alow,aupp));

            ipanel++;
        }
    }

    TCanvas *c3 = new TCanvas("c3", "", 900, 900);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.1);
    for (int it=0; it<nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];
        for (int ia=0; ia<nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            if (it==0 && ia==0)
                hX2copy[it][ia]->Draw("P");
            else
                hX2copy[it][ia]->Draw("P SAME");
        }
    }
    leg2->Draw("SAME");
}

//******************************************************************************
//******************************************************************************

void LoadData()
{
    fin = TFile::Open(input.Data(), "READ");
    for (int it=0; it<nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];
        for (int ia=0; ia<nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;

            hX1[it][ia] = (TH1D*)fin->Get(Form("xFractions/hX1[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hX2[it][ia] = (TH1D*)fin->Get(Form("xFractions/hX2[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hX1X2[it][ia] = (TH2D*)fin->Get(Form("xFractions/hX1X2[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

            // Config histos
            hX1[it][ia]->SetTitle(";Fraction x;1/N dx/dN");
            hX1[it][ia]->SetMarkerStyle(kFullCircle);
            hX1[it][ia]->SetMarkerSize(0.8);
            hX1[it][ia]->SetMarkerColor(kBlack);
            hX1[it][ia]->SetLineColor(kBlack);
            hX1[it][ia]->SetLineWidth(2);
            hX1[it][ia]->GetXaxis()->SetLabelSize(0.062);
            hX1[it][ia]->GetYaxis()->SetLabelSize(0.062);
            hX1[it][ia]->GetXaxis()->SetTitleSize(0.072);
            hX1[it][ia]->GetYaxis()->SetTitleSize(0.072);
            hX1[it][ia]->GetYaxis()->SetMaxDigits(3);
            hX1[it][ia]->Scale(1./hX1[it][ia]->GetEntries());
            hX1[it][ia]->GetXaxis()->SetRangeUser(2e-7, 1.);
            hX1[it][ia]->GetYaxis()->SetRangeUser(-0.0001, 0.051);
            hX1[it][ia]->GetYaxis()->SetTitleOffset(0.8);
            hX1[it][ia]->GetXaxis()->SetTitleOffset(0.91);

            hX2[it][ia]->SetTitle(";Fraction x;1/N dx/dN");
            hX2[it][ia]->SetMarkerStyle(kFullCircle);
            hX2[it][ia]->SetMarkerSize(0.8);
            hX2[it][ia]->SetMarkerColor(kRed);
            hX2[it][ia]->SetLineColor(kRed);
            hX2[it][ia]->SetLineWidth(2);
            hX2[it][ia]->GetXaxis()->SetLabelSize(0.062);
            hX2[it][ia]->GetYaxis()->SetLabelSize(0.062);
            hX2[it][ia]->GetXaxis()->SetTitleSize(0.072);
            hX2[it][ia]->GetYaxis()->SetTitleSize(0.072);
            hX2[it][ia]->GetYaxis()->SetMaxDigits(3);
            hX2[it][ia]->Scale(1./hX2[it][ia]->GetEntries());

            hX1X2[it][ia]->SetTitle(";x_{1};x_{2}");
            hX1X2[it][ia]->Rebin2D(2,2);
            hX1X2[it][ia]->GetYaxis()->SetMaxDigits(2);
            hX1X2[it][ia]->GetXaxis()->SetLabelSize(0.062);
            hX1X2[it][ia]->GetYaxis()->SetLabelSize(0.062);
            hX1X2[it][ia]->GetXaxis()->SetTitleSize(0.072);
            hX1X2[it][ia]->GetYaxis()->SetTitleSize(0.072);
            hX1X2[it][ia]->Scale(1./hX1X2[it][ia]->GetEntries());
            hX1X2[it][ia]->GetXaxis()->SetRangeUser(2e-7, 1.);
            hX1X2[it][ia]->GetYaxis()->SetRangeUser(2e-7, 1.);
            hX1X2[it][ia]->GetYaxis()->SetTitleOffset(0.8);
            hX1X2[it][ia]->GetXaxis()->SetTitleOffset(0.91);


            hX2copy[it][ia] = (TH1D*)hX2[it][ia]->Clone(Form("hX2copy[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hX2copy[it][ia]->Rebin(4);
            hX2copy[it][ia]->GetYaxis()->SetRangeUser(0.000001, 1.);
            hX2copy[it][ia]->GetXaxis()->SetLabelSize(0.032);
            hX2copy[it][ia]->GetYaxis()->SetLabelSize(0.032);
            hX2copy[it][ia]->GetXaxis()->SetTitleSize(0.042);
            hX2copy[it][ia]->GetYaxis()->SetTitleSize(0.042);
            hX2copy[it][ia]->SetMarkerSize(1.);
            hX2copy[it][ia]->SetMarkerStyle(mStyle[it][ia]);
            hX2copy[it][ia]->SetMarkerColor(mColor[it][ia]);
            hX2copy[it][ia]->SetLineColor(mColor[it][ia]);
            hX2copy[it][ia]->GetYaxis()->SetTitleOffset(1.1);
            hX2copy[it][ia]->GetXaxis()->SetTitleOffset(0.91);
        }
    }
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
    gStyle->SetLineWidth(2);
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

void redrawBorder()
{
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
   l.DrawLine(2e-7, gPad->GetUymax(), 1., gPad->GetUymax());
   l.DrawLine(1., gPad->GetUymin(), 1., gPad->GetUymax());
   l.DrawLine(2e-7, gPad->GetUymin(), 1., gPad->GetUymin());
}