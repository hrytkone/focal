void LoadData();
void ConfigPlots();
void DrawData();
void SetStyle(Bool_t graypalette);

const int nTriggBins = 3;
const int nAssocBins = 4;
const double triggPt[nTriggBins+1] = {1.0, 2.0, 2.5, 3.0};
const double assocPt[nAssocBins+1] = {0.5, 1.0, 1.5, 2.0, 2.5};

const int nset = 2;
const int mColor[nset] = {2,1};
const TString infiles[nset] = {
    "/home/heimarry/Simulations/focal/STAR-data/output_pAu.root",
    "/home/heimarry/Simulations/focal/STAR-data/output_pp.root"
};


TFile *fin[nset];
TH2D *hCorr[nset][nTriggBins][nAssocBins];
TH1D *hTrigg[nset];
TH1D *hProjX[nset][nTriggBins][nAssocBins];
TCanvas *cCorr[nTriggBins][nAssocBins];

int nRealTrigg[nset][nTriggBins] = {0};
double maxval[nTriggBins][nAssocBins] = {0};
double minval[nTriggBins][nAssocBins] = {0};

void CompareRealPi0()
{
    SetStyle(0);
	LoadData();
    ConfigPlots();
	DrawData();
}

void LoadData()
{
	for (int iset = 0; iset < nset; iset++) {
        fin[iset] = TFile::Open(infiles[iset]);
    }

	for (int iset = 0; iset < nset; iset++) {
        hTrigg[iset] = (TH1D*)fin[iset]->Get(Form("hRealTriggCounter"));
		for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            nRealTrigg[iset][itrigg] = hTrigg[iset]->GetBinContent(itrigg+1);
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
			for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;

				hCorr[iset][itrigg][iassoc] = (TH2D*)fin[iset]->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hProjX[iset][itrigg][iassoc] = hCorr[iset][itrigg][iassoc]->ProjectionX(Form("ProjX[%d][%4.1f,%4.1f][%4.1f,%4.1f]",iset,tlow,tupp,alow,aupp),0,-1,"e");

                maxval[itrigg][iassoc] = -1.;
                minval[itrigg][iassoc] = 10000000000.;
			}
		}
	}
}

void ConfigPlots()
{
    for (int iset = 0; iset < nset; iset++) {
        for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
            double tlow = triggPt[itrigg];
            double tupp = triggPt[itrigg+1];
            for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
                double alow = assocPt[iassoc];
                double aupp = assocPt[iassoc+1];

                if (tlow < aupp) continue;

                hCorr[iset][itrigg][iassoc] = (TH2D*)fin[iset]->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
                hProjX[iset][itrigg][iassoc] = hCorr[iset][itrigg][iassoc]->ProjectionX(Form("ProjX[%d][%4.1f,%4.1f][%4.1f,%4.1f]",iset,tlow,tupp,alow,aupp),0,-1,"e");
                hProjX[iset][itrigg][iassoc]->SetMarkerStyle(8);
                hProjX[iset][itrigg][iassoc]->SetLineColor(mColor[iset]);
                hProjX[iset][itrigg][iassoc]->SetMarkerColor(mColor[iset]);
                hProjX[iset][itrigg][iassoc]->GetYaxis()->SetMaxDigits(2);
                hProjX[iset][itrigg][iassoc]->Rebin(4);
                hProjX[iset][itrigg][iassoc]->Scale(1./nRealTrigg[iset][itrigg], "width");

                double maximum = hProjX[iset][itrigg][iassoc]->GetBinContent(hProjX[iset][itrigg][iassoc]->GetMaximumBin());
                if (maximum>maxval[itrigg][iassoc])
                    maxval[itrigg][iassoc] = maximum;
                double minimum = hProjX[iset][itrigg][iassoc]->GetBinContent(hProjX[iset][itrigg][iassoc]->GetMinimumBin());
                if (minimum<minval[itrigg][iassoc])
                    minval[itrigg][iassoc] = minimum;
            }
        }
    }
}

void DrawData()
{
    cout << "Number of triggers : " << endl;
    //cCorr[0][0] = new TCanvas(Form("cCorr_%d_%d", 0, 0), Form("cCorr_%d_%d", 0, 0), 600, 600);
    //cCorr[0][0]->Divide(2,1);
    //cCorr[0][0]->cd(1);
    //hProjX[0][0][0]->Draw("");
    //cCorr[0][0]->cd(2);
    //hProjX[1][0][0]->Draw("");

    TLegend *leg = new TLegend(0.68, 0.75, 0.88, 0.85);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
    leg->AddEntry(hProjX[0][0][0], "pAu", "pe");
    leg->AddEntry(hProjX[1][0][0], "pp", "pe");

    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
        for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

            cCorr[itrigg][iassoc] = new TCanvas(Form("cCorr_%d_%d", itrigg, iassoc), Form("cCorr_%d_%d", itrigg, iassoc), 600, 600);

            for (int iset = 0; iset < nset; iset++) {
                cout << Form("\tset %d, [%4.1f,%4.1f][%4.1f,%4.1f] : ntrigg=",iset,tlow,tupp,alow,aupp) << nRealTrigg[iset][itrigg] << "  npair=" << hCorr[iset][itrigg][iassoc]->GetEntries() << endl;
                //hProjX[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(TMath::Pi()/2., 2*TMath::Pi());
                hProjX[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(minval[itrigg][iassoc]-0.1*minval[itrigg][iassoc], maxval[itrigg][iassoc]+0.1*maxval[itrigg][iassoc]);
                hProjX[iset][itrigg][iassoc]->Draw("SAME");
            }
            leg->Draw("SAME");
            cCorr[itrigg][iassoc]->SaveAs(Form("corr_%d_%d.pdf", itrigg, iassoc));
        }
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
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
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
