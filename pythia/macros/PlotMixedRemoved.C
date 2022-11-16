const int nTriggBins = 4;
const int nAssocBins = 5;
const double triggPt[nTriggBins+1] = {4.0, 8.0, 10.0, 15.0, 20.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0, 8.0, 10.0, 15.0};

int itrigg = 0;
int iassoc = 0;
double edgecut = 0.8;

int binx = 0, biny = 0;
double alpha = 0;

TFile *fin;
TH2D *hCorrMassMass[nTriggBins][nAssocBins];
TH2D *hCorrSideMass[nTriggBins][nAssocBins];
TH2D *hCorrMassSide[nTriggBins][nAssocBins];
TH2D *hCorrSideSide[nTriggBins][nAssocBins];

TH2D *hCorrMassMassMixed[nTriggBins][nAssocBins];
TH2D *hCorrMassSideMixed[nTriggBins][nAssocBins];
TH2D *hCorrSideMassMixed[nTriggBins][nAssocBins];
TH2D *hCorrSideSideMixed[nTriggBins][nAssocBins];

TH2D *hCorrMassMassDiv[nTriggBins][nAssocBins];
TH2D *hCorrMassSideDiv[nTriggBins][nAssocBins];
TH2D *hCorrSideMassDiv[nTriggBins][nAssocBins];
TH2D *hCorrSideSideDiv[nTriggBins][nAssocBins];

TH2D *hCorrMixedSum;

TCanvas *canvasMM, *canvasMS, *canvasSM, *canvasSS;

void LoadData(TString input);

void PlotMixedRemoved(TString input="input.root")
{
    gStyle->SetOptStat(0);

    LoadData(input);

    hCorrMassMass[itrigg][iassoc]->SetTitle(Form("Same event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrMassSide[itrigg][iassoc]->SetTitle(Form("Same event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrSideMass[itrigg][iassoc]->SetTitle(Form("Same event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrSideSide[itrigg][iassoc]->SetTitle(Form("Same event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));

    hCorrMassMassMixed[itrigg][iassoc]->SetTitle(Form("Mixed event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrMassSideMixed[itrigg][iassoc]->SetTitle(Form("Mixed event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrSideMassMixed[itrigg][iassoc]->SetTitle(Form("Mixed event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrSideSideMixed[itrigg][iassoc]->SetTitle(Form("Mixed event [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));

    hCorrMassMassDiv[itrigg][iassoc]->SetTitle(Form("Corrected [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrMassSideDiv[itrigg][iassoc]->SetTitle(Form("Corrected [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrSideMassDiv[itrigg][iassoc]->SetTitle(Form("Corrected [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));
    hCorrSideSideDiv[itrigg][iassoc]->SetTitle(Form("Corrected [%.1f,%.1f][%.1f,%.1f];#Delta#phi;#Delta#eta", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));

    // Mass-mass
    canvasMM = new TCanvas("canvasMM", "canvasMM", 1800, 600);
    canvasMM->Divide(3,1);
    canvasMM->cd(1);
    gPad->SetLogz();
    hCorrMassMass[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrMassMass[itrigg][iassoc]->Draw("LEGO2");
    canvasMM->cd(2);
    //gPad->SetLogz();
    binx = hCorrMassMassMixed[itrigg][iassoc]->GetXaxis()->FindBin(0.);
    biny = hCorrMassMassMixed[itrigg][iassoc]->GetYaxis()->FindBin(0.);
    alpha = 1./hCorrMassMassMixed[itrigg][iassoc]->GetBinContent(binx, biny);
    hCorrMassMassMixed[itrigg][iassoc]->Scale(alpha);
    hCorrMassMassMixed[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrMassMassMixed[itrigg][iassoc]->Draw("LEGO2");
    canvasMM->cd(3);
    gPad->SetLogz();
    hCorrMassMassDiv[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrMassMassDiv[itrigg][iassoc]->Divide(hCorrMassMassMixed[itrigg][iassoc]);
    hCorrMassMassDiv[itrigg][iassoc]->Draw("LEGO2");
    canvasMM->SaveAs(Form("mass-mass[%.1f,%.1f][%.1f,%.1f].png", triggPt[itrigg], triggPt[itrigg+1], assocPt[iassoc], assocPt[iassoc+1]));

    // Mass-side
    canvasMS = new TCanvas("canvasMS", "canvasMS", 1800, 600);
    canvasMS->Divide(3,1);
    canvasMS->cd(1);
    gPad->SetLogz();
    hCorrMassSide[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrMassSide[itrigg][iassoc]->Draw("LEGO2");
    canvasMS->cd(2);
    //gPad->SetLogz();
    binx = hCorrMassSideMixed[itrigg][iassoc]->GetXaxis()->FindBin(0.);
    biny = hCorrMassSideMixed[itrigg][iassoc]->GetYaxis()->FindBin(0.);
    alpha = 1./hCorrMassSideMixed[itrigg][iassoc]->GetBinContent(binx, biny);
    hCorrMassSideMixed[itrigg][iassoc]->Scale(alpha);
    hCorrMassSideMixed[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrMassSideMixed[itrigg][iassoc]->Draw("LEGO2");
    canvasMS->cd(3);
    gPad->SetLogz();
    hCorrMassSideDiv[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrMassSideDiv[itrigg][iassoc]->Divide(hCorrMassSideMixed[itrigg][iassoc]);
    hCorrMassSideDiv[itrigg][iassoc]->Draw("LEGO2");
    canvasMS->SaveAs(Form("mass-side[%.1f,%.1f][%.1f,%.1f].png", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));

    // Side-mass
    canvasSM = new TCanvas("canvasSM", "canvasSM", 1800, 600);
    canvasSM->Divide(3,1);
    canvasSM->cd(1);
    gPad->SetLogz();
    hCorrSideMass[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrSideMass[itrigg][iassoc]->Draw("LEGO2");
    canvasSM->cd(2);
    //gPad->SetLogz();
    binx = hCorrSideMassMixed[itrigg][iassoc]->GetXaxis()->FindBin(0.);
    biny = hCorrSideMassMixed[itrigg][iassoc]->GetYaxis()->FindBin(0.);
    alpha = 1./hCorrSideMassMixed[itrigg][iassoc]->GetBinContent(binx, biny);
    hCorrSideMassMixed[itrigg][iassoc]->Scale(alpha);
    hCorrSideMassMixed[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrSideMassMixed[itrigg][iassoc]->Draw("LEGO2");
    canvasSM->cd(3);
    gPad->SetLogz();
    hCorrSideMassDiv[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrSideMassDiv[itrigg][iassoc]->Divide(hCorrSideMassMixed[itrigg][iassoc]);
    hCorrSideMassDiv[itrigg][iassoc]->Draw("LEGO2");
    canvasSM->SaveAs(Form("side-mass[%.1f,%.1f][%.1f,%.1f].png", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));

    // Side-side
    canvasSS = new TCanvas("canvasSS", "canvasSS", 1800, 600);
    canvasSS->Divide(3,1);
    canvasSS->cd(1);
    gPad->SetLogz();
    hCorrSideSide[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrSideSide[itrigg][iassoc]->Draw("LEGO2");
    canvasSS->cd(2);
    //gPad->SetLogz();
    binx = hCorrSideSideMixed[itrigg][iassoc]->GetXaxis()->FindBin(0.);
    biny = hCorrSideSideMixed[itrigg][iassoc]->GetYaxis()->FindBin(0.);
    alpha = 1./hCorrSideSideMixed[itrigg][iassoc]->GetBinContent(binx, biny);
    hCorrSideSideMixed[itrigg][iassoc]->Scale(alpha);
    hCorrSideSideMixed[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrSideSideMixed[itrigg][iassoc]->Draw("LEGO2");
    canvasSS->cd(3);
    gPad->SetLogz();
    hCorrSideSideDiv[itrigg][iassoc]->GetYaxis()->SetRangeUser(-edgecut, edgecut);
    hCorrSideSideDiv[itrigg][iassoc]->Divide(hCorrSideSideMixed[itrigg][iassoc]);
    hCorrSideSideDiv[itrigg][iassoc]->Draw("LEGO2");
    canvasSS->SaveAs(Form("side-side[%.1f,%.1f][%.1f,%.1f].png", triggPt[itrigg], triggPt[itrigg+1], assocPt[itrigg], assocPt[itrigg+1]));

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

            hCorrMassMass[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassMass[itrigg][iassoc]->Rebin2D(6,1);
            hCorrMassMassDiv[itrigg][iassoc] = (TH2D*)hCorrMassMass[itrigg][iassoc]->Clone(Form("hCorrMassMassDiv_%d_%d", itrigg, iassoc));

            hCorrMassSide[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassSide[itrigg][iassoc]->Rebin2D(6,1);
            hCorrMassSideDiv[itrigg][iassoc] = (TH2D*)hCorrMassSide[itrigg][iassoc]->Clone(Form("hCorrMassSideDiv_%d_%d", itrigg, iassoc));

            hCorrSideMass[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideMass[itrigg][iassoc]->Rebin2D(6,1);
            hCorrSideMassDiv[itrigg][iassoc] = (TH2D*)hCorrSideMass[itrigg][iassoc]->Clone(Form("hCorrSideMassDiv_%d_%d", itrigg, iassoc));

            hCorrSideSide[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideSide[itrigg][iassoc]->Rebin2D(6,1);
            hCorrSideSideDiv[itrigg][iassoc] = (TH2D*)hCorrSideSide[itrigg][iassoc]->Clone(Form("hCorrSideSideDiv_%d_%d", itrigg, iassoc));

			hCorrMassMassMixed[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrMassMass/hCorrMassMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassMassMixed[itrigg][iassoc]->Rebin2D(6,1);

            hCorrMassSideMixed[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrMassSide/hCorrMassSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrMassSideMixed[itrigg][iassoc]->Rebin2D(6,1);

            hCorrSideMassMixed[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrSideMass/hCorrSideMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideMassMixed[itrigg][iassoc]->Rebin2D(6,1);

            hCorrSideSideMixed[itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrSideSide/hCorrSideSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            hCorrSideSideMixed[itrigg][iassoc]->Rebin2D(6,1);

            if (itrigg==0 && iassoc==0) {
                hCorrMixedSum = (TH2D*)hCorrMassMassMixed[itrigg][iassoc]->Clone("hCorrMixedSum");
            } else {
                hCorrMixedSum->Add(hCorrMassMassMixed[itrigg][iassoc]);
            }
		}
	}
}
