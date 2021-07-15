#include "include/Filipad.h"
#include "include/rootcommon.h"

void LoadData();
void Draw();

double lowx=-20.;
double highx=20.;
double ly = 1e-5;
double hy = 100.;
double lowIAA = -0.2;
double highIAA = 5;

const int nTriggBins = 4;
const int nAssocBins = 4;
double triggPt[nTriggBins+1] = {1.0, 2.0, 4.0, 8.0, 20.0};
double assocPt[nAssocBins+1] = {0.5, 1.0, 2.0, 3.0, 4.0};

const int nset = 1;
TString infiles[nset] = {
	"analyse_test_1m_focal.root"
};
TFile *fin[nset];
TH1D *hCorr[nset][nTriggBins][nAssocBins];

void ComparisonPlot()
{
	LoadData();
	DrawFiliPad();
}

void LoadData()
{
	for (int iset = 0; iset < nset; iset++)
		fin[iset] = TFile::Open(infiles[iset]);

	for (int iset = 0; iset < nset; iset++) {
		for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
			for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
				hCorr[iset][itrigg][iassoc] = (TH1D*)fin[iset]->Get(Form("hCorr%d:%d",itrigg,iassoc));
			}		
		}
	}
}

void DrawFiliPad()
{
	int padID = 1;
	Filipad *fpad = new Filipad(padID, 1.1, 0.5, 100, 100, 0.7, 1);
	fpad->Draw();

	// Upper pad
	TPad *p = fpad->GetPad(1);
	p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
	hset(*hCorr[0][2][2], "#delta#phi", "Counts", 1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
	hCorr[0][2][2]->Draw();
}