const int nTriggBins = 3;
const int nAssocBins = 4;
const double triggPt[nTriggBins+1] = {2.0, 4.0, 8.0, 20.0};
const double assocPt[nAssocBins+1] = {1.0, 1.5, 2.0, 3.0, 4.0};

TFile *fin;
TH2D *hCorrMixed[nTriggBins][nAssocBins];
TCanvas *canvas[nTriggBins][nAssocBins];
TCanvas *canvas_px[nTriggBins][nAssocBins];

void LoadData(TString input);

void PlotMixedEvent(TString input="input.root")
{
    //gStyle->SetOptStat(0);

    LoadData(input);
    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
		for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

			canvas[itrigg][iassoc] = new TCanvas(Form("c_%d_%d", itrigg, iassoc), "", 600, 600);
            hCorrMixed[itrigg][iassoc]->Draw("LEGO2Z");

			canvas_px[itrigg][iassoc] = new TCanvas(Form("cpx_%d_%d", itrigg, iassoc), "", 600, 600);
            hCorrMixed[itrigg][iassoc]->ProjectionX()->Draw("P");

		}
	}
}

//******************************************************************************
//******************************************************************************

void LoadData(TString input)
{
    fin = TFile::Open(input.Data());
	for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
		for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

			hCorrMixed[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrMassMass/hCorrMassMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMixed[itrigg][iassoc]->Rebin2D(8);
		}
	}
}
