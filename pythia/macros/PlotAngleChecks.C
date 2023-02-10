
//TString infile = "/home/heimarry/Simulations/focal/analysis_output/2023-01-10_pp-focal_oangle-divided_asym-1.root";
//TString infile = "/home/heimarry/Simulations/focal/analysis_output/2023-01-09_pp-focal_oangle-divided_asym-08.root";
//TString infile = "/home/heimarry/Simulations/focal/analysis_output/2023-01-10_pp-focal_oangle-divided_asym-05.root";
TString infile = "/home/heimarry/Simulations/focal/analysis_output_2/2023-02-09_pp-focal_smaller-mw_no-asym.root";


const int nAssocBins = 3;
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0, 8.0};

TFile *fin;

TH2D* hMassAsymTrue[nAssocBins];
TH2D* hMassAsymDecay[nAssocBins];
TH2D* hMassAsymMix[nAssocBins];
TH2D* hMassAsymNotDecay[nAssocBins];

TH2D* hMassOpeningAngleTrue[nAssocBins];
TH2D* hMassOpeningAngleDecay[nAssocBins];
TH2D* hMassOpeningAngleMix[nAssocBins];
TH2D* hMassOpeningAngleNotDecay[nAssocBins];
TH2D* hMassOpeningAngleSum[nAssocBins];

TCanvas *cAsym[nAssocBins];
TCanvas *cOpeningAngle[nAssocBins];

void LoadData();
void ConfigHistos();
void PlotAsym();
void PlotOpeningAngle();

void PlotAngleChecks()
{
    gStyle->SetOptStat(0);

    LoadData();
    ConfigHistos();
    PlotAsym();
    PlotOpeningAngle();
}

//******************************************************************************
//******************************************************************************

void LoadData()
{
	fin = TFile::Open(infile);
	for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
        double alow = assocPt[iassoc];
        double aupp = assocPt[iassoc+1];

        hMassAsymTrue[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassAsymTrue[%4.1f,%4.1f]",alow,aupp));
        hMassAsymDecay[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassAsymDecay[%4.1f,%4.1f]",alow,aupp));
        hMassAsymMix[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassAsymMix[%4.1f,%4.1f]",alow,aupp));
        hMassAsymNotDecay[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassAsymNotDecay[%4.1f,%4.1f]",alow,aupp));

        hMassOpeningAngleTrue[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassOpeningAngleTrue[%4.1f,%4.1f]",alow,aupp));
        hMassOpeningAngleDecay[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassOpeningAngleDecay[%4.1f,%4.1f]",alow,aupp));
        hMassOpeningAngleMix[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassOpeningAngleMix[%4.1f,%4.1f]",alow,aupp));
        hMassOpeningAngleNotDecay[iassoc]  = (TH2D*)fin->Get(Form("Asymmetry/hMassOpeningAngleNotDecay[%4.1f,%4.1f]",alow,aupp));
	}
}

void ConfigHistos()
{
    for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
        double alow = assocPt[iassoc];
        double aupp = assocPt[iassoc+1];

        hMassOpeningAngleSum[iassoc] = (TH2D*)hMassOpeningAngleTrue[iassoc]->Clone(Form("hMassOpeningAngleSum[%4.1f,%4.1f]",alow,aupp));
        hMassOpeningAngleSum[iassoc]->Add(hMassOpeningAngleDecay[iassoc]);
        hMassOpeningAngleSum[iassoc]->Add(hMassOpeningAngleMix[iassoc]);
        hMassOpeningAngleSum[iassoc]->Add(hMassOpeningAngleNotDecay[iassoc]);

        hMassOpeningAngleTrue[iassoc]->Divide(hMassOpeningAngleSum[iassoc]);
        hMassOpeningAngleDecay[iassoc]->Divide(hMassOpeningAngleSum[iassoc]);
        hMassOpeningAngleMix[iassoc]->Divide(hMassOpeningAngleSum[iassoc]);
        hMassOpeningAngleNotDecay[iassoc]->Divide(hMassOpeningAngleSum[iassoc]);

        hMassOpeningAngleTrue[iassoc]->GetYaxis()->SetRangeUser(0., 0.05);
        hMassOpeningAngleDecay[iassoc]->GetYaxis()->SetRangeUser(0., 0.05);
        hMassOpeningAngleMix[iassoc]->GetYaxis()->SetRangeUser(0., 0.05);
        hMassOpeningAngleNotDecay[iassoc]->GetYaxis()->SetRangeUser(0., 0.05);
    }
}

void PlotAsym()
{
	for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
        cAsym[iassoc] = new TCanvas(Form("cAsym_%d", iassoc), "", 800, 600);
        cAsym[iassoc]->Divide(2,2);
        cAsym[iassoc]->cd(1);
        hMassAsymTrue[iassoc]->Draw("COLZ");
        cAsym[iassoc]->cd(2);
        hMassAsymDecay[iassoc]->Draw("COLZ");
        cAsym[iassoc]->cd(3);
        hMassAsymMix[iassoc]->Draw("COLZ");
        cAsym[iassoc]->cd(4);
        hMassAsymNotDecay[iassoc]->Draw("COLZ");
    }
}

void PlotOpeningAngle()
{
    for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
        cOpeningAngle[iassoc] = new TCanvas(Form("cOpeningAngle_%d", iassoc), "", 600, 600);
        cOpeningAngle[iassoc]->Divide(2,2);
        cOpeningAngle[iassoc]->cd(1);
        hMassOpeningAngleTrue[iassoc]->Draw("COL");
        cOpeningAngle[iassoc]->cd(2);
        hMassOpeningAngleDecay[iassoc]->Draw("COL");
        cOpeningAngle[iassoc]->cd(3);
        hMassOpeningAngleMix[iassoc]->Draw("COL");
        cOpeningAngle[iassoc]->cd(4);
        hMassOpeningAngleNotDecay[iassoc]->Draw("COL");
    }
}
