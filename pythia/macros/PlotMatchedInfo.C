void LoadData(TString input);
void ConfigHistos();
void DrawHistos();

const int nhisto = 3;
const TString hNames[nhisto] = {"hPtMatched", "hEtaMatched", "hPhiMatched"};

TFile *fin;
TH2D *histo[nhisto];
TCanvas *canvas[nhisto];

void PlotMatchedInfo(TString input)
{
    gStyle->SetOptStat(0);

	LoadData(input);
    //ConfigHistos();
	DrawHistos();
}

//************************************************************************************************
//************************************************************************************************

void LoadData(TString input)
{
	fin = TFile::Open(input);
    for (int i = 0; i < nhisto; i++)
        histo[i]  = (TH2D*)fin->Get(hNames[i].Data());
}

void ConfigHistos()
{

}

void DrawHistos()
{
    for (int i = 0; i < nhisto; i++) {
        canvas[i] = new TCanvas(Form("canvas_%s", hNames[i].Data()), "", 600, 600);
        histo[i]->Draw("COLZ");
    }
}
