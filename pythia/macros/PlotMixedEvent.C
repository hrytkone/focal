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

TH2D *hCorrMixedSum;
TH1D *hCorrMixedSumX;
TH1D *hCorrMixedSumY;

TCanvas *canvas[nTriggBins][nAssocBins];
TCanvas *canvas_px_py[nTriggBins][nAssocBins];

TCanvas *canvasSum;
TCanvas *canvasSum_px_py;

void LoadData(TString input);

void PlotMixedEvent(TString input="input.root")
{
    gStyle->SetOptStat(0);

    LoadData(input);

    canvasSum = new TCanvas("csum", "", 600, 600);
    hCorrMixedSum->Draw("LEGO2Z");

    canvasSum_px_py = new TCanvas("canvasSum_px_py", "canvasSum_px_py", 1200, 600);
    canvasSum_px_py->Divide(2,1);

    canvasSum_px_py->cd(1);
    hCorrMixedSumX = hCorrMixedSum->ProjectionX();
    int maxBin = hCorrMixedSumX->GetMaximumBin();
    double rangeMax = hCorrMixedSumX->GetBinContent(maxBin) + hCorrMixedSumX->GetBinContent(maxBin);
    hCorrMixedSumX->GetYaxis()->SetRangeUser(0., rangeMax);
    hCorrMixedSumX->SetMarkerStyle(kOpenCircle);
    hCorrMixedSumX->SetMarkerColor(kBlack);
    hCorrMixedSumX->SetLineColor(kBlack);
    hCorrMixedSumX->SetTitle(Form("%s;#phi;counts", hCorrMixedSum->GetName()));
    hCorrMixedSumX->Draw("P");

    canvasSum_px_py->cd(2);
    hCorrMixedSumY = hCorrMixedSum->ProjectionY();
    maxBin = hCorrMixedSumY->GetMaximumBin();
    rangeMax = hCorrMixedSumY->GetBinContent(maxBin) + 0.7*hCorrMixedSumY->GetBinContent(maxBin);
    hCorrMixedSumY->GetYaxis()->SetRangeUser(0., rangeMax);
    hCorrMixedSumY->SetMarkerStyle(kOpenCircle);
    hCorrMixedSumY->SetMarkerColor(kBlack);
    hCorrMixedSumY->SetLineColor(kBlack);
    hCorrMixedSumY->SetTitle(Form("%s;#phi;counts", hCorrMixedSum->GetName()));
    hCorrMixedSumY->Draw("P");

    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
		for (int iassoc = 0; iassoc < nTriggBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

			canvas[itrigg][iassoc] = new TCanvas(Form("c_%d_%d", itrigg, iassoc), "", 600, 600);
            hCorrMixed[itrigg][iassoc]->Draw("LEGO2Z");
            //canvas[itrigg][iassoc]->SaveAs(Form("corr2d_500k_2.png");

			canvas_px_py[itrigg][iassoc] = new TCanvas(Form("cpx_%d_%d", itrigg, iassoc), "", 1200, 600);
            canvas_px_py[itrigg][iassoc]->Divide(2,1);

            canvas_px_py[itrigg][iassoc]->cd(1);
            hCorrMixedProjX[itrigg][iassoc] = hCorrMixed[itrigg][iassoc]->ProjectionX();
            maxBin = hCorrMixedProjX[itrigg][iassoc]->GetMaximumBin();
            rangeMax = hCorrMixedProjX[itrigg][iassoc]->GetBinContent(maxBin) + hCorrMixedProjX[itrigg][iassoc]->GetBinContent(maxBin);
            hCorrMixedProjX[itrigg][iassoc]->GetYaxis()->SetRangeUser(0., rangeMax);
            hCorrMixedProjX[itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
            hCorrMixedProjX[itrigg][iassoc]->SetMarkerColor(kBlack);
            hCorrMixedProjX[itrigg][iassoc]->SetLineColor(kBlack);
            hCorrMixedProjX[itrigg][iassoc]->SetTitle(Form("[%.1f,%.1f][%.1f,%.1f];#Delta#phi;counts", triggPt[itrigg], triggPt[itrigg+1], assocPt[iassoc], assocPt[iassoc+1]));
            hCorrMixedProjX[itrigg][iassoc]->Draw("P");

            canvas_px_py[itrigg][iassoc]->cd(2);
            hCorrMixedProjY[itrigg][iassoc] = hCorrMixed[itrigg][iassoc]->ProjectionY();
            maxBin = hCorrMixedProjY[itrigg][iassoc]->GetMaximumBin();
            rangeMax = hCorrMixedProjY[itrigg][iassoc]->GetBinContent(maxBin) + 0.7*hCorrMixedProjY[itrigg][iassoc]->GetBinContent(maxBin);
            hCorrMixedProjY[itrigg][iassoc]->GetYaxis()->SetRangeUser(0., rangeMax);
            hCorrMixedProjY[itrigg][iassoc]->SetMarkerStyle(kOpenCircle);
            hCorrMixedProjY[itrigg][iassoc]->SetMarkerColor(kBlack);
            hCorrMixedProjY[itrigg][iassoc]->SetLineColor(kBlack);
            hCorrMixedProjY[itrigg][iassoc]->SetTitle(";#Delta#eta;");
            hCorrMixedProjY[itrigg][iassoc]->Draw("P");
            canvas_px_py[itrigg][iassoc]->SaveAs(Form("proj_[%.1f,%.1f][%.1f,%.1f].png", triggPt[itrigg], triggPt[itrigg+1], assocPt[iassoc], assocPt[iassoc+1]));
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

            if (itrigg==0 && iassoc==0) {
                hCorrMixedSum = (TH2D*)hCorrMixed[itrigg][iassoc]->Clone("hCorrMixedSum");
            } else {
                hCorrMixedSum->Add(hCorrMixed[itrigg][iassoc]);
            }
            cout << "bin width" << endl;
            cout << "\tdelta phi : " << hCorrMixed[itrigg][iassoc]->GetXaxis()->GetBinWidth(0) << endl;
            cout << "\tdelta eta : " << hCorrMixed[itrigg][iassoc]->GetYaxis()->GetBinWidth(0) << endl;
		}
	}
}
