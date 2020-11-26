const int nTriggBins = 8;
double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};
//const int nTriggBins = 4;
//double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 15.0, 20.0};

const int nAssocBins = 7;
double  assocPt[nAssocBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};
//const int nAssocBins = 6;
//double  assocPt[nAssocBins+1] = {1.0, 2.0, 5.0, 6.0, 8.0, 10.0, 15.0};

void PlotCorrelations(TString sInputName = "output.root", int triggbin = 0, int assocbin = 0)
{
    TFile *fIn = TFile::Open(sInputName);

    TH2D *hCorrMid[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrFor[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrChargedMid[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrChargedFor[nTriggBins+1][nAssocBins+1];

    for (int it=0; it<nTriggBins; it++) {
        for (int ia=0; ia<nAssocBins; ia++) {
            hCorrMid[it][ia] = (TH2D*)fIn->Get(Form("hCorrMid%d:%d", it, ia));
            hCorrFor[it][ia] = (TH2D*)fIn->Get(Form("hCorrFor%d:%d", it, ia));
            hCorrChargedMid[it][ia] = (TH2D*)fIn->Get(Form("hCorrChargedMid%d:%d", it, ia));
            hCorrChargedFor[it][ia] = (TH2D*)fIn->Get(Form("hCorrChargedFor%d:%d", it, ia));
        }
    }

    gStyle->SetOptStat(0);

    double nTriggMid = hCorrMid[triggbin][assocbin]->GetEntries();
    double nTriggFor = hCorrFor[triggbin][assocbin]->GetEntries();
    TH1D *hCorrMidProj = hCorrMid[triggbin][assocbin]->ProjectionX();
    TH1D *hCorrForProj = hCorrFor[triggbin][assocbin]->ProjectionX();
    hCorrMidProj->Scale(1.0/nTriggMid, "width");
    hCorrForProj->Scale(1.0/nTriggFor, "width");

    hCorrMidProj->SetTitle("; #Delta#phi; 1/N_{trigg }dN/d#phi");
    hCorrMidProj->SetMarkerStyle(20);
    hCorrMidProj->SetMarkerColor(2);
    hCorrMidProj->GetYaxis()->SetRangeUser(0.01, 2.0);
    
    hCorrForProj->SetMarkerStyle(20);
    hCorrForProj->SetMarkerColor(4);
    
    TLegend *leg = new TLegend(0.68, 0.68, 0.88, 0.85);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetHeader(Form("#splitline{%.1f < p_{T,t} < %.1f GeV/c}{%.1f < p_{T,a} < %.1f GeV/c}", triggPt[triggbin], triggPt[triggbin+1], assocPt[assocbin], assocPt[assocbin+1]));
    leg->AddEntry(hCorrMidProj, "-1.0 < #eta < 1.0", "p");
    leg->AddEntry(hCorrForProj, "3.2 < #eta < 5.8", "p");
    
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->SetLogy(1);
    hCorrMidProj->Draw("P");
    hCorrForProj->Draw("P SAME");
    leg->Draw("SAME");
}
