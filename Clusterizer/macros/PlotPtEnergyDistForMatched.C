TFile *fin;
TH2D *hEnergyPhoton;
TH2D *hEnergyPhotonFromDecay;
TH2D *hEnergyChargedPion;
TH2D *hEnergyElectron;
TH1D *hEnergyPhotonCluster;
TH1D *hEnergyPhotonClusterFromDecay;
TH1D *hEnergyChargedPionCluster;
TH1D *hEnergyElectronCluster;

TH2D *hPtPhoton;
TH2D *hPtPhotonFromDecay;
TH2D *hPtChargedPion;
TH2D *hPtElectron;
TH1D *hPtPhotonCluster;
TH1D *hPtPhotonClusterFromDecay;
TH1D *hPtChargedPionCluster;
TH1D *hPtElectronCluster;

TH1D *hEnergyHCALPhoton;
TH1D *hEnergyHCALPhotonFromDecay;
TH1D *hEnergyHCALChargedPion;
TH1D *hEnergyHCALElectron;

TCanvas *cEnergyCluster;
TCanvas *cPtCluster;
TCanvas *cEnergyHCAL;

void LoadData(TString input);
void ConfigPlots();
void Plot();

void PlotPtEnergyDistForMatched(TString input)
{
    LoadData(input);
    ConfigPlots();
    Plot();
}

//______________________________________________________________________________
void LoadData(TString input)
{
    fin = TFile::Open(input.Data());
    hEnergyPhoton = (TH2D*)fin->Get("hEnergyPhoton");
    hEnergyPhotonFromDecay = (TH2D*)fin->Get("hEnergyPhotonFromDecay");
    hEnergyChargedPion = (TH2D*)fin->Get("hEnergyChargedPion");
    hEnergyElectron = (TH2D*)fin->Get("hEnergyElectron");
    hEnergyPhotonCluster = hEnergyPhoton->ProjectionY();
    hEnergyPhotonClusterFromDecay = hEnergyPhotonFromDecay->ProjectionY();
    hEnergyChargedPionCluster = hEnergyChargedPion->ProjectionY();
    hEnergyElectronCluster = hEnergyElectron->ProjectionY();

    hPtPhoton = (TH2D*)fin->Get("hPtPhoton");
    hPtPhotonFromDecay = (TH2D*)fin->Get("hPtPhotonFromDecay");
    hPtChargedPion = (TH2D*)fin->Get("hPtChargedPion");
    hPtElectron = (TH2D*)fin->Get("hPtElectron");
    hPtPhotonCluster = hPtPhoton->ProjectionY();
    hPtPhotonClusterFromDecay = hPtPhotonFromDecay->ProjectionY();
    hPtChargedPionCluster = hPtChargedPion->ProjectionY();
    hPtElectronCluster = hPtElectron->ProjectionY();

    hEnergyHCALPhoton = (TH1D*)fin->Get("hEnergyHCALPhoton");
    hEnergyHCALPhotonFromDecay = (TH1D*)fin->Get("hEnergyHCALPhotonFromDecay");
    hEnergyHCALChargedPion = (TH1D*)fin->Get("hEnergyHCALChargedPion");
    hEnergyHCALElectron = (TH1D*)fin->Get("hEnergyHCALElectron");
}

void ConfigPlots()
{
    hEnergyPhotonCluster->SetMarkerColor(kBlack);
    hEnergyPhotonCluster->SetMarkerStyle(kCircle);
    hEnergyPhotonClusterFromDecay->SetMarkerColor(kGray);
    hEnergyPhotonClusterFromDecay->SetMarkerStyle(kCircle);
    hEnergyElectronCluster->SetMarkerColor(kBlue);
    hEnergyElectronCluster->SetMarkerStyle(kCircle);
    hEnergyChargedPionCluster->SetMarkerColor(kRed);
    hEnergyChargedPionCluster->SetMarkerStyle(kCircle);

    hPtPhotonCluster->SetMarkerColor(kBlack);
    hPtPhotonCluster->SetMarkerStyle(kCircle);
    hPtPhotonClusterFromDecay->SetMarkerColor(kGray);
    hPtPhotonClusterFromDecay->SetMarkerStyle(kCircle);
    hPtElectronCluster->SetMarkerColor(kBlue);
    hPtElectronCluster->SetMarkerStyle(kCircle);
    hPtChargedPionCluster->SetMarkerColor(kRed);
    hPtChargedPionCluster->SetMarkerStyle(kCircle);

    hEnergyHCALPhoton->SetMarkerColor(kBlack);
    hEnergyHCALPhoton->SetMarkerStyle(kCircle);
    hEnergyHCALPhotonFromDecay->SetMarkerColor(kGray);
    hEnergyHCALPhotonFromDecay->SetMarkerStyle(kCircle);
    hEnergyHCALElectron->SetMarkerColor(kBlue);
    hEnergyHCALElectron->SetMarkerStyle(kCircle);
    hEnergyHCALChargedPion->SetMarkerColor(kRed);
    hEnergyHCALChargedPion->SetMarkerStyle(kCircle);
}

void Plot()
{
    cEnergyCluster = new TCanvas("c1", "", 600, 600);
    cEnergyCluster->SetLogy();
    hEnergyPhotonCluster->Draw("P");
    hEnergyPhotonClusterFromDecay->Draw("P SAME");
    hEnergyElectronCluster->Draw("P SAME");
    hEnergyChargedPionCluster->Draw("P SAME");

    cPtCluster = new TCanvas("c2", "", 600, 600);
    cPtCluster->SetLogy();
    hPtPhotonCluster->Draw("P");
    hPtPhotonClusterFromDecay->Draw("P SAME");
    hPtElectronCluster->Draw("P SAME");
    hPtChargedPionCluster->Draw("P SAME");

    cEnergyHCAL = new TCanvas("c3", "", 600, 600);
    cEnergyHCAL->SetLogy();
    hEnergyHCALPhoton->Draw("P");
    hEnergyHCALPhotonFromDecay->Draw("P SAME");
    hEnergyHCALElectron->Draw("P SAME");
    hEnergyHCALChargedPion->Draw("P SAME");
}
