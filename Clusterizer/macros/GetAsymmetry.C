const int nset = 3;
const TString filename[nset] = {
    "/home/heimarry/Simulations/focal/analysis_output_2/2023-01-27_pp-focal_asym-08.root",
    "/home/heimarry/Simulations/focal/analysis_output_2/2023-01-31_pp-focal_asym-05.root",
    "/home/heimarry/Simulations/focal/analysis_output_2/2023-01-31_gAnalysis_asym-08_etacut-02.root"
};

const TString histname[nset] = {
    "hEnergyAsymTrue",
    "hEnergyAsymRec",
    "hEnergyAsymRec"
};

const TString houtname[nset] = {
    "hAsymTrue",
    "hAsymRecPythia",
    "hAsymRecGeant"
};

TFile *fin[nset];
TH1D *hAsymmetry[nset];
TH1D *hCounter[nset];

void LoadData()
{
    for (int i = 0; i < nset; i++) {
        fin[i] = TFile::Open(filename[i].Data());
        hCounter[i] = (TH1D*)fin[i]->Get("hCounter");
        hAsymmetry[i] = (TH1D*)fin[i]->Get(Form("%s", histname[i].Data()));
        hAsymmetry[i]->Scale(1./hCounter[i]->GetBinContent(1));
    }
}

void GetAsymmetry()
{
    LoadData();
    TFile *fout = TFile::Open("asymmetry_output.root", "RECREATE");
    for (int i = 0; i < nset; i++) {
        TString newname = Form("%s", houtname[i].Data());
        fout->WriteObject(hAsymmetry[i], newname.Data());
    }
    fout->Close();
}
