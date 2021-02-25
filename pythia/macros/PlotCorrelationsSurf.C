const int nTriggBins = 8;
double  triggPt[nTriggBins+1] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};

const int nAssocBins = 7;
double  assocPt[nAssocBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0};

void PlotCorrelationsSurf(TString sInputName = "output.root")
{

    int mintrigg = 2, maxtrigg = 3;
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
    for (int triggbin=mintrigg; triggbin<maxtrigg; triggbin++) { 
        for (int assocbin=minassoc; assocbin<maxassoc; assocbin++) {

            c1->cd(icanvas);
            icanvas++;

            double nTriggMid = hCorrMid[triggbin][assocbin]->GetEntries();
            double nTriggFor = hCorrFor[triggbin][assocbin]->GetEntries();
            double nTriggMidRec = hCorrPion0RecMid[triggbin][assocbin]->GetEntries();
            double nTriggForRec = hCorrPion0RecFor[triggbin][assocbin]->GetEntries();
    
            hCorrFor[triggbin][assocbin]->Draw("SURF1");
        }
    }
}
