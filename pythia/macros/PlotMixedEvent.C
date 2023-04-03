const int nset = 4;
const int nTriggBins = 2;
const int nAssocBins = 2;
const double triggPt[nTriggBins+1] = {4.0, 8.0, 16.0};
const double assocPt[nAssocBins+1] = {2.0, 3.0, 4.0};

TFile *fin;
TH2D *hCorrMixed[nset][nTriggBins][nAssocBins];
TH1D *hCorrMixedProjX[nset][nTriggBins][nAssocBins];
TH1D *hCorrMixedProjY[nset][nTriggBins][nAssocBins];

TCanvas *canvas[nTriggBins][nAssocBins];
TCanvas *canvas2[nTriggBins][nAssocBins];

const TString label[nset] = {
    Form("mass-mass"),
    Form("mass-side"),
    Form("side-mass"),
    Form("side-side")
};

void LoadData(TString input);
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin);
double XtoPad(double x);
double YtoPad(double y);
void SetStyle(Bool_t graypalette);

void PlotMixedEvent(TString input="input.root")
{
    SetStyle(0);
    LoadData(input);

    for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
		for (int iassoc = 0; iassoc < nTriggBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

			canvas[itrigg][iassoc] = new TCanvas(Form("c_%d_%d", itrigg, iassoc), "", 600, 600);
			canvas2[itrigg][iassoc] = new TCanvas(Form("c2_%d_%d", itrigg, iassoc), "", 600, 600);

            // Number of PADS
            const Int_t Nx = 2;
            const Int_t Ny = 2;

            // Margins
            Float_t lMargin = 0.11;
            Float_t rMargin = 0.02;
            Float_t bMargin = 0.11;
            Float_t tMargin = 0.02;

            // Canvas setup
            CanvasPartition(canvas[itrigg][iassoc],Nx,Ny,lMargin,rMargin,bMargin,tMargin);
            CanvasPartition(canvas2[itrigg][iassoc],Nx,Ny,lMargin,rMargin,bMargin,tMargin);

            int labelsize = 16;
            TPad *pad[Nx][Ny];
            int iset = 0;
            int comp[nset] = {2, 0, 3, 1};
            for (Int_t i = 0; i < Nx; i++) {
                for (Int_t j = 0; j < Ny; j++) {
                    canvas[itrigg][iassoc]->cd(0);

                    // Get the pads previously created.
                    pad[i][j] = (TPad*) canvas[itrigg][iassoc]->FindObject(TString::Format("pad_%d_%d",i,j).Data());
                    pad[i][j]->Draw();
                    pad[i][j]->SetFillStyle(4000);
                    pad[i][j]->SetFrameFillStyle(4000);
                    pad[i][j]->cd();

                    hCorrMixedProjX[comp[iset]][itrigg][iassoc]->Draw("HIST E");
                    if (i==1 && j==1) {
                        TLatex l;
                        l.SetTextFont(43);
                        l.SetTextSize(labelsize);
                        l.DrawLatexNDC(XtoPad(0.4), YtoPad(0.85), Form("p_{T,t} = [%0.1f, %0.1f] GeV/c", triggPt[itrigg], triggPt[itrigg+1]));
                        l.DrawLatexNDC(XtoPad(0.4), YtoPad(0.78), Form("p_{T,a} = [%0.1f, %0.1f] GeV/c", assocPt[iassoc], assocPt[iassoc+1]));
                    }

                    TLatex lbl;
                    lbl.SetTextFont(43);
                    lbl.SetTextSize(labelsize);
                    lbl.DrawLatexNDC(XtoPad(0.1), YtoPad(0.85), label[comp[iset]]);

                    iset++;
                }
            }

            iset = 0;
            for (Int_t i = 0; i < Nx; i++) {
                for (Int_t j = 0; j < Ny; j++) {
                    canvas2[itrigg][iassoc]->cd(0);

                    // Get the pads previously created.
                    pad[i][j] = (TPad*) canvas2[itrigg][iassoc]->FindObject(TString::Format("pad_%d_%d",i,j).Data());
                    pad[i][j]->Draw();
                    pad[i][j]->SetFillStyle(4000);
                    pad[i][j]->SetFrameFillStyle(4000);
                    pad[i][j]->cd();

                    hCorrMixedProjY[comp[iset]][itrigg][iassoc]->Draw("HIST E");
                    if (i==1 && j==1) {
                        TLatex l;
                        l.SetTextFont(43);
                        l.SetTextSize(labelsize);
                        l.DrawLatexNDC(XtoPad(0.4), YtoPad(0.88), Form("p_{T,t} = [%0.1f, %0.1f] GeV/c", triggPt[itrigg], triggPt[itrigg+1]));
                        l.DrawLatexNDC(XtoPad(0.4), YtoPad(0.81), Form("p_{T,a} = [%0.1f, %0.1f] GeV/c", assocPt[iassoc], assocPt[iassoc+1]));
                    }

                    TLatex lbl;
                    lbl.SetTextFont(43);
                    lbl.SetTextSize(labelsize);
                    lbl.DrawLatexNDC(XtoPad(0.1), YtoPad(0.88), label[comp[iset]]);

                    iset++;
                }
            }
		}
	}
}

//******************************************************************************
//******************************************************************************

void LoadData(TString input)
{
    fin = TFile::Open(input.Data());
	for (int itrigg = 0; itrigg < nTriggBins; itrigg++) {
        double tlow = triggPt[itrigg];
        double tupp = triggPt[itrigg+1];
		for (int iassoc = 0; iassoc < nAssocBins; iassoc++) {
            double alow = assocPt[iassoc];
            double aupp = assocPt[iassoc+1];

            if (tlow < aupp) continue;

			hCorrMixed[0][itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrMassMass/hCorrMassMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
			hCorrMixed[1][itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrMassSide/hCorrMassSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
			hCorrMixed[2][itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrSideMass/hCorrSideMassMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
			hCorrMixed[3][itrigg][iassoc] = (TH2D*)fin->Get(Form("CorrSideSide/hCorrSideSideMixed[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));

            //hCorrMixed[0][itrigg][iassoc]->Rebin2D(4,2);
            //hCorrMixed[1][itrigg][iassoc]->Rebin2D(4,2);
            //hCorrMixed[2][itrigg][iassoc]->Rebin2D(4,2);
            //hCorrMixed[3][itrigg][iassoc]->Rebin2D(4,2);

            int rebin = 3;
            for (int iset = 0; iset < nset; iset++) {
                hCorrMixedProjX[iset][itrigg][iassoc] = hCorrMixed[iset][itrigg][iassoc]->ProjectionX();
                hCorrMixedProjX[iset][itrigg][iassoc]->Rebin(rebin);
                hCorrMixedProjX[iset][itrigg][iassoc]->Scale(1./hCorrMixedProjX[iset][itrigg][iassoc]->GetEntries(), "width");
                hCorrMixedProjX[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.024, 0.29);
                hCorrMixedProjX[iset][itrigg][iassoc]->SetLineColor(kBlack);
                hCorrMixedProjX[iset][itrigg][iassoc]->SetLineWidth(2);
                hCorrMixedProjX[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrMixedProjX[iset][itrigg][iassoc]->SetMarkerStyle(kFullCircle);
                hCorrMixedProjX[iset][itrigg][iassoc]->SetMarkerSize(0);


                hCorrMixedProjX[iset][itrigg][iassoc]->GetXaxis()->SetTitleFont(43);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetXaxis()->SetLabelFont(43);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetXaxis()->SetTitleSize(24);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetXaxis()->SetLabelSize(20);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetXaxis()->SetTitleOffset(1.5);

                hCorrMixedProjX[iset][itrigg][iassoc]->GetYaxis()->SetTitleFont(43);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetYaxis()->SetLabelFont(43);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetYaxis()->SetTitleSize(24);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetYaxis()->SetLabelSize(20);
                hCorrMixedProjX[iset][itrigg][iassoc]->GetYaxis()->SetTitleOffset(2.5);

                if (iset == 0) hCorrMixedProjX[iset][itrigg][iassoc]->SetTitle(";;1/N dN/d(#Delta#phi)");
                else if (iset == nset-1) hCorrMixedProjX[iset][itrigg][iassoc]->SetTitle(";#Delta#phi;");
                else hCorrMixedProjX[iset][itrigg][iassoc]->SetTitle(";;");

                hCorrMixedProjY[iset][itrigg][iassoc] = hCorrMixed[iset][itrigg][iassoc]->ProjectionY();
                hCorrMixedProjY[iset][itrigg][iassoc]->Rebin(rebin);
                hCorrMixedProjY[iset][itrigg][iassoc]->Scale(1./hCorrMixedProjY[iset][itrigg][iassoc]->GetEntries(), "width");
                hCorrMixedProjY[iset][itrigg][iassoc]->GetXaxis()->SetRangeUser(-2.49, 2.49);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetYaxis()->SetRangeUser(-0.05, 0.84);
                hCorrMixedProjY[iset][itrigg][iassoc]->SetLineColor(kBlack);
                hCorrMixedProjY[iset][itrigg][iassoc]->SetLineWidth(2);
                hCorrMixedProjY[iset][itrigg][iassoc]->SetMarkerColor(kBlack);
                hCorrMixedProjY[iset][itrigg][iassoc]->SetMarkerStyle(kFullCircle);
                hCorrMixedProjY[iset][itrigg][iassoc]->SetMarkerSize(0);


                hCorrMixedProjY[iset][itrigg][iassoc]->GetXaxis()->SetTitleFont(43);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetXaxis()->SetLabelFont(43);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetXaxis()->SetTitleSize(24);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetXaxis()->SetLabelSize(20);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetXaxis()->SetTitleOffset(1.5);

                hCorrMixedProjY[iset][itrigg][iassoc]->GetYaxis()->SetTitleFont(43);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetYaxis()->SetLabelFont(43);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetYaxis()->SetTitleSize(24);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetYaxis()->SetLabelSize(20);
                hCorrMixedProjY[iset][itrigg][iassoc]->GetYaxis()->SetTitleOffset(2.5);

                if (iset == 0) hCorrMixedProjY[iset][itrigg][iassoc]->SetTitle(";;1/N dN/d(#Delta#eta)");
                else if (iset == nset-1) hCorrMixedProjY[iset][itrigg][iassoc]->SetTitle(";#Delta#eta;");
                else hCorrMixedProjY[iset][itrigg][iassoc]->SetTitle(";;");
            }
		}
	}
}

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   for (Int_t i=0;i<Nx;i++) {

      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j=0;j<Ny;j++) {

         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         auto name = TString::Format("pad_%d_%d",i,j);
         auto pad = (TPad*) C->FindObject(name.Data());
         if (pad) delete pad;
         pad = new TPad(name.Data(),"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);

         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);

         pad->Draw();
      }
   }
}

double XtoPad(double x)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double pw = xu-xl;
   double lm = gPad->GetLeftMargin();
   double rm = gPad->GetRightMargin();
   double fw = pw-pw*lm-pw*rm;
   return (x*fw+pw*lm)/pw;
}

double YtoPad(double y)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double ph = yu-yl;
   double tm = gPad->GetTopMargin();
   double bm = gPad->GetBottomMargin();
   double fh = ph-ph*bm-ph*tm;
   return (y*fh+bm*ph)/ph;
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;
    gStyle->SetTextFont(43);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
}
