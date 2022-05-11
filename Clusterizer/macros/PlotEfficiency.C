const TString filename = "efficiency_small-bins_no-asym-cut.root";

TFile *fin;

TH2D *hEtaPtTrue;
TH2D *hEtaPtRec;
TH2D *hEfficiency;
TH1D *hEfficiency_py;
TH1D *hEtaPtTrue_py;
TH1D *hEtaPtRec_py;

TH1D *hEPhotonTrue;
TH1D *hEPhotonRec;

void LoadData();

void PlotEfficiency()
{
    gStyle->SetOptStat(0);

    LoadData();
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    c1->SetLogy();
    c1->SetLogz();
    hEfficiency->SetTitle(";#eta;p_{T}");
    hEfficiency->Draw("COLZ");

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    c2->SetLogx();
    hEfficiency_py->SetTitle(";p_{T};efficiency");
    hEfficiency_py->Draw("PE");

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
    c3->SetLogx();
    hEtaPtTrue_py->SetTitle(";p_{T};efficiency");
    hEtaPtTrue_py->Draw("PE");
    hEtaPtRec_py->Draw("PE SAME");

    TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
    hEPhotonTrue->Draw("PE");
    hEPhotonRec->Draw("PE SAME");
}

void LoadData()
{
    fin = TFile::Open(filename.Data());
    hEtaPtTrue = (TH2D*)fin->Get("hEtaETrue");
    hEtaPtRec = (TH2D*)fin->Get("hEtaERec");
    hEfficiency = (TH2D*)hEtaPtRec->Clone("hEfficiency");
    hEfficiency->Divide(hEtaPtTrue);
    hEtaPtTrue_py = hEtaPtTrue->ProjectionY();
    hEtaPtTrue_py->SetLineColor(kBlack);
    hEtaPtRec_py = hEtaPtRec->ProjectionY();
    hEtaPtRec_py->SetLineColor(kRed);
    hEfficiency_py = (TH1D*)hEtaPtRec_py->Clone("hEfficiency");
    hEfficiency_py->Divide(hEtaPtTrue_py);
    hEfficiency_py->Rebin(10);

    hEPhotonTrue = (TH1D*)fin->Get("hEPhotonTrue");
    hEPhotonRec = (TH1D*)fin->Get("hEPhotonCluster");
    hEPhotonTrue->SetLineColor(kRed);
    hEPhotonRec->SetLineColor(kBlack);
}
