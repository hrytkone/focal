const int nset = 6;
//const TString filenames[nset] = {
//    "matching_pythia_cut-5.root", // pT
//    "matching_pythia_cut-3.root", // pT, alpha
//    "matching_pythia_cut-6.root", // pT, alpha, Eclust
//    "matching_pythia_cut-1.root", // pT, eta, Eclust
//    "matching_pythia_cut-0.root", // pT, eta, alpha
//    "matching_pythia_cut-2.root"  // pT, eta, alpha, Eclust
//};

//const TString filenames[nset] = {
//    "matching_pythia_asymcut-0.root",
//    "matching_pythia_asymcut-1.root",
//    "matching_pythia_asymcut-2.root",
//    "matching_pythia_asymcut-3.root",
//    "matching_pythia_asymcut-4.root",
//    "matching_pythia_asymcut-5.root"
//};

//const TString filenames[nset] = {
//    "matching_pythia_pt-0.root",
//    "matching_pythia_pt-1.root",
//    "matching_pythia_pt-2.root",
//    "matching_pythia_pt-3.root",
//    "matching_pythia_pt-4.root",
//    "matching_pythia_pt-5.root"
//};

//const TString filenames[nset] = {
//    "matching_pythia_asymcut-pt-0.root",
//    "matching_pythia_asymcut-pt-1.root",
//    "matching_pythia_asymcut-pt-2.root",
//    "matching_pythia_asymcut-pt-3.root",
//    "matching_pythia_asymcut-pt-4.root",
//    "matching_pythia_asymcut-pt-5.root"
//};

const TString filenames[nset] = {
    "matching_pythia_asymcut-etacut-pt-0.root",
    "matching_pythia_asymcut-etacut-pt-1.root",
    "matching_pythia_asymcut-etacut-pt-2.root",
    "matching_pythia_asymcut-etacut-pt-3.root",
    "matching_pythia_asymcut-etacut-pt-4.root",
    "matching_pythia_asymcut-etacut-pt-5.root"
};

//------------------------------------------------------------------------------

//const TString legHeader[nset] = {
//    "5 GeV < p_{T} < 12 GeV",
//    "5 GeV < p_{T} < 12 GeV",
//    "5 GeV < p_{T} < 12 GeV",
//    "5 GeV < p_{T} < 12 GeV",
//    "5 GeV < p_{T} < 12 GeV",
//    "5 GeV < p_{T} < 12 GeV"
//};

const TString legHeader[nset] = {
    "1 GeV < p_{T} < 2 GeV",
    "2 GeV < p_{T} < 3 GeV",
    "3 GeV < p_{T} < 5 GeV",
    "5 GeV < p_{T} < 8 GeV",
    "8 GeV < p_{T} < 12 GeV",
    "12 GeV < p_{T} < 20 GeV"
};

//------------------------------------------------------------------------------

const TString legEntry1[nset] = {
    "#alpha < 0.5",
    "#alpha < 0.5",
    "#alpha < 0.5",
    "#alpha < 0.5",
    "#alpha < 0.5",
    "#alpha < 0.5"
};

//const TString legEntry1[nset] = {
//    "#alpha < 0.4",
//    "#alpha < 0.5",
//    "#alpha < 0.6",
//    "#alpha < 0.7",
//    "#alpha < 0.8",
//    "#alpha < 0.9"
//};

//const TString legEntry1[nset] = {
//    " ",
//    "#alpha < 0.5",
//    "#alpha < 0.5",
//    "4 < #eta < 5",
//    "4 < #eta < 5",
//    "4 < #eta < 5"
//};

//------------------------------------------------------------------------------

const TString legEntry2[nset] = {
    "4 < #eta < 5",
    "4 < #eta < 5",
    "4 < #eta < 5",
    "4 < #eta < 5",
    "4 < #eta < 5",
    "4 < #eta < 5"
};

//const TString legEntry2[nset] = {
//    " ",
//    " ",
//    "E_{clust} > 5 GeV",
//    "E_{clust} > 5 GeV",
//    "#alpha < 0.5",
//    "#alpha < 0.5"
//};

//------------------------------------------------------------------------------

//const TString legEntry3[nset] = {
//    " ",
//    " ",
//    " ",
//    " ",
//    " ",
//    "E_{clust} > 5 GeV"
//};

//------------------------------------------------------------------------------

//const double xi[nset] = {0.45, 0.45, 0.45, 0.45, 0.45, 0.45};
//const double xf[nset] = {0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
//const double yi[nset] = {0.25, 0.55, 0.55, 0.55, 0.55, 0.55};
//const double yf[nset] = {0.45, 0.85, 0.85, 0.85, 0.85, 0.85};

const double xi[nset] = {0.45, 0.45, 0.45, 0.45, 0.45, 0.45};
const double xf[nset] = {0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
const double yi[nset] = {0.20, 0.65, 0.65, 0.65, 0.65, 0.65};
const double yf[nset] = {0.40, 0.85, 0.85, 0.85, 0.85, 0.85};

//------------------------------------------------------------------------------

TFile *fin[nset];
TH1D *hMassCluster[nset];
TLegend *leg[nset];

void LoadData();
void CreateLegends();

void PlotMassPeaks()
{
    gStyle->SetOptStat(0);

    LoadData();
    CreateLegends();
    TCanvas *c1 = new TCanvas("c1", "c1", 900, 600);
    c1->Divide(3,2);
    for (int i=0; i<nset; i++) {
        c1->cd(i+1);
        gPad->SetLeftMargin(0.1);
        gPad->SetBottomMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        hMassCluster[i]->SetTitle(";M_{#gamma#gamma};counts");
        hMassCluster[i]->Draw("PE");
        leg[i]->Draw("SAME");
    }
}

void LoadData()
{
    for (int ifile=0; ifile<nset; ifile++) {
        fin[ifile] = TFile::Open(filenames[ifile].Data());
        hMassCluster[ifile] = (TH1D*)fin[ifile]->Get("hMassCluster");
        hMassCluster[ifile]->Rebin();
        hMassCluster[ifile]->SetMarkerStyle(kDot);
        hMassCluster[ifile]->SetMarkerSize(.5);
        hMassCluster[ifile]->SetMarkerColor(kBlack);
        hMassCluster[ifile]->SetLineColor(kBlack);
        hMassCluster[ifile]->GetYaxis()->SetMaxDigits(3);
    }
}

void CreateLegends()
{
    for (int i=0; i<nset; i++) {
        leg[i] = new TLegend(xi[i], yi[i], xf[i], yf[i]);
        leg[i]->SetFillStyle(0); leg[i]->SetBorderSize(0); leg[i]->SetTextSize(gStyle->GetTextSize()*1.);
        leg[i]->SetHeader(legHeader[i].Data());
        leg[i]->AddEntry(hMassCluster[i], legEntry1[i], "");
        leg[i]->AddEntry(hMassCluster[i], legEntry2[i], "");
        //leg[i]->AddEntry(hMassCluster[i], legEntry3[i], "");
    }
}
