const int nTriggBins = 8;
double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};
//const int nTriggBins = 4;
//double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 15.0, 20.0};

const int nAssocBins = 7;
double  assocPt[nAssocBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};
//const int nAssocBins = 6;
//double  assocPt[nAssocBins+1] = {1.0, 2.0, 5.0, 6.0, 8.0, 10.0, 15.0};

void PlotCorrelations(TString sInputName = "output.root")
{

    int mintrigg = 1, maxtrigg = 2;
    int minassoc = 0, maxassoc = 6;
    const double mSize = .5;

    TFile *fIn = TFile::Open(sInputName);

    TH2D *hCorrMid[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrFor[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrChargedMid[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrChargedFor[nTriggBins+1][nAssocBins+1];
    
    TH2D *hCorrPion0RecMid[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrPion0RecFor[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrPion0SideMid[nTriggBins+1][nAssocBins+1];
    TH2D *hCorrPion0SideFor[nTriggBins+1][nAssocBins+1];

    for (int it=0; it<nTriggBins; it++) {
        for (int ia=0; ia<nAssocBins; ia++) {
            hCorrMid[it][ia] = (TH2D*)fIn->Get(Form("hCorrMid%d:%d", it, ia));
            hCorrFor[it][ia] = (TH2D*)fIn->Get(Form("hCorrFor%d:%d", it, ia));
            hCorrChargedMid[it][ia] = (TH2D*)fIn->Get(Form("hCorrChargedMid%d:%d", it, ia));
            hCorrChargedFor[it][ia] = (TH2D*)fIn->Get(Form("hCorrChargedFor%d:%d", it, ia));
            hCorrPion0RecMid[it][ia] = (TH2D*)fIn->Get(Form("hCorrPion0RecMid%d:%d", it, ia));
            hCorrPion0RecFor[it][ia] = (TH2D*)fIn->Get(Form("hCorrPion0RecFor%d:%d", it, ia));
            hCorrPion0SideMid[it][ia] = (TH2D*)fIn->Get(Form("hCorrPion0SideMid%d:%d", it, ia));
            hCorrPion0SideFor[it][ia] = (TH2D*)fIn->Get(Form("hCorrPion0SideFor%d:%d", it, ia));
        }
    }

    gStyle->SetOptStat(0);
    
    TLegend *leg = new TLegend(0.38, 0.68, 0.58, 0.85);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
    
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
    //c1->Divide(maxtrigg-mintrigg, maxassoc-minassoc);
    c1->Divide(3, 2);
   
    int icanvas = 1;
    bool bPlotLeg = true;
    for (int triggbin=mintrigg; triggbin<maxtrigg; triggbin++) { 
    for (int assocbin=minassoc; assocbin<maxassoc; assocbin++) {

        c1->cd(icanvas);
        icanvas++;

        double nTriggMid = hCorrMid[triggbin][assocbin]->GetEntries();
        double nTriggFor = hCorrFor[triggbin][assocbin]->GetEntries();
        double nTriggMidRec = hCorrPion0RecMid[triggbin][assocbin]->GetEntries();
        double nTriggForRec = hCorrPion0RecFor[triggbin][assocbin]->GetEntries();
    
        TH1D *hCorrMidProj = hCorrMid[triggbin][assocbin]->ProjectionX();
        TH1D *hCorrForProj = hCorrFor[triggbin][assocbin]->ProjectionX();
        hCorrMidProj->Scale(1.0/nTriggMid, "width");
        hCorrForProj->Scale(1.0/nTriggFor, "width");

        TH1D *hCorrRecMidProj = hCorrPion0RecMid[triggbin][assocbin]->ProjectionX();
        TH1D *hCorrRecForProj = hCorrPion0RecFor[triggbin][assocbin]->ProjectionX();
        hCorrRecMidProj->Scale(1.0/nTriggMidRec, "width");
        hCorrRecForProj->Scale(1.0/nTriggForRec, "width");
    
        hCorrMidProj->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; 1/N_{trigg}dN/d#phi", triggPt[triggbin], triggPt[triggbin+1], assocPt[assocbin], assocPt[assocbin+1]));
        hCorrMidProj->SetMarkerStyle(20);
        hCorrMidProj->SetMarkerColor(2);
        hCorrMidProj->SetMarkerSize(mSize);
        hCorrMidProj->GetYaxis()->SetRangeUser(0.0, 2.1);
    
        hCorrForProj->SetMarkerStyle(20);
        hCorrForProj->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; 1/N_{trigg}dN/d#phi", triggPt[triggbin], triggPt[triggbin+1], assocPt[assocbin], assocPt[assocbin+1]));
        hCorrForProj->SetMarkerColor(2);
        hCorrForProj->SetMarkerSize(mSize);
 
        hCorrRecMidProj->SetMarkerStyle(20);
        hCorrRecMidProj->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; 1/N_{trigg}dN/d#phi", triggPt[triggbin], triggPt[triggbin+1], assocPt[assocbin], assocPt[assocbin+1]));
        hCorrRecMidProj->SetMarkerColor(2);
        hCorrRecMidProj->SetMarkerSize(mSize);
 
        hCorrRecForProj->SetMarkerStyle(20);
        hCorrRecForProj->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; 1/N_{trigg}dN/d#phi", triggPt[triggbin], triggPt[triggbin+1], assocPt[assocbin], assocPt[assocbin+1]));
        hCorrRecForProj->SetMarkerColor(4);
        hCorrRecForProj->SetMarkerSize(mSize);

        hCorrForProj->Draw("P SAME");
        hCorrMidProj->Draw("P");
        //hCorrRecForProj->Draw("P");
        //hCorrRecMidProj->Draw("P SAME");
        if (bPlotLeg) {
            leg->AddEntry(hCorrMidProj, "PYTHIA #pi^{0}, -1.0 < #eta < 1.0", "p");
            leg->AddEntry(hCorrForProj, "PYTHIA #pi^{0}, 3.2 < #eta < 5.8", "p");
            //leg->AddEntry(hCorrRecMidProj, "rec. #pi^{0}, -1.0 < #eta < 1.0", "p");
            //leg->AddEntry(hCorrRecForProj, "rec. #pi^{0}, 3.2 < #eta < 5.8", "p");
            leg->Draw("SAME");
            bPlotLeg = false;
        }
    }
    }
}
