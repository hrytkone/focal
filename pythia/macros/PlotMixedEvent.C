const int nTriggBins = 4;
const int nAssocBins = 5;
//const double triggPt[nTriggBins+1] = {4.0, 8.0, 20.0};
//const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};
const double triggPt[nTriggBins+1] = {4.0, 8.0, 10.0, 15.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0, 8.0, 10.0, 15.0};

TFile *fin;
TH2D *hCorrMixed[nTriggBins][nAssocBins];
TH1D *hCorrMixedProjX[nTriggBins][nAssocBins];
TH1D *hCorrMixedProjY[nTriggBins][nAssocBins];
TCanvas *canvas[nTriggBins][nAssocBins];
TCanvas *canvas_px_py[nTriggBins][nAssocBins];

void LoadData(TString input);

void PlotMixedEvent(TString input="input.root")
{
    //gStyle->SetOptStat(0);

    LoadData(input);
    for (int itrigg = 0; itrigg < 1; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
		for (int iassoc = 0; iassoc < 1; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

			canvas[itrigg][iassoc] = new TCanvas(Form("c_%d_%d", itrigg, iassoc), "", 600, 600);
            hCorrMixed[itrigg][iassoc]->Draw("LEGO2Z");

			canvas_px_py[itrigg][iassoc] = new TCanvas(Form("cpx_%d_%d", itrigg, iassoc), "", 1200, 600);
            canvas_px_py[itrigg][iassoc]->Divide(2,1);

            canvas_px_py[itrigg][iassoc]->cd(1);
            hCorrMixedProjX[itrigg][iassoc] = hCorrMixed[itrigg][iassoc]->ProjectionX();
            int maxBin = hCorrMixedProjX[itrigg][iassoc]->GetMaximumBin();
            double rangeMax = hCorrMixedProjX[itrigg][iassoc]->GetBinContent(maxBin) + hCorrMixedProjX[itrigg][iassoc]->GetBinContent(maxBin);
            hCorrMixedProjX[itrigg][iassoc]->GetYaxis()->SetRangeUser(0., rangeMax);
            hCorrMixedProjX[itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
            hCorrMixedProjX[itrigg][iassoc]->SetMarkerColor(kBlack);
            hCorrMixedProjX[itrigg][iassoc]->SetLineColor(kBlack);
            hCorrMixedProjX[itrigg][iassoc]->SetTitle(";#phi;counts");
            hCorrMixedProjX[itrigg][iassoc]->Draw("P");

            canvas_px_py[itrigg][iassoc]->cd(2);
            hCorrMixedProjY[itrigg][iassoc] = hCorrMixed[itrigg][iassoc]->ProjectionY();
            maxBin = hCorrMixedProjY[itrigg][iassoc]->GetMaximumBin();
            rangeMax = hCorrMixedProjY[itrigg][iassoc]->GetBinContent(maxBin) + 0.7*hCorrMixedProjY[itrigg][iassoc]->GetBinContent(maxBin);
            hCorrMixedProjY[itrigg][iassoc]->GetYaxis()->SetRangeUser(0., rangeMax);
            hCorrMixedProjY[itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
            hCorrMixedProjY[itrigg][iassoc]->SetMarkerColor(kBlack);
            hCorrMixedProjY[itrigg][iassoc]->SetLineColor(kBlack);
            hCorrMixedProjY[itrigg][iassoc]->SetTitle(";#eta;");
            hCorrMixedProjY[itrigg][iassoc]->Draw("P");
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
            hCorrMixed[itrigg][iassoc]->Rebin2D(4,2);
            cout << "bin width" << endl;
            cout << "\tdelta phi : " << hCorrMixed[itrigg][iassoc]->GetXaxis()->GetBinWidth(0) << endl;
            cout << "\tdelta eta : " << hCorrMixed[itrigg][iassoc]->GetYaxis()->GetBinWidth(0) << endl;
		}
	}
}
