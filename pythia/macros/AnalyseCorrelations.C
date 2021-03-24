const int nTriggBins = 4;
double  triggPt[nTriggBins+1] = {1.0, 2.0, 4.0, 8.0, 20.0};

const int nAssocBins = 4;
double  assocPt[nAssocBins+1] = {0.5, 1.0, 2.0, 3.0, 4.0};

double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);

void AnalyseCorrelations(TString sInputName = "output.root")
{

    double brGammaCh = 0.988;

    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    gStyle->SetStatY(0.85);                
    gStyle->SetStatX(0.85);                
    gStyle->SetStatW(0.2);                
    gStyle->SetStatH(0.1);                
    
    TFile *fIn = TFile::Open(sInputName);   

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hCounter = (TH1D*)fIn->Get("hCounter");
    int nEvent = hCounter->GetBinContent(1);
    std::cout << "Input file contains " << nEvent << " events, proceed to analyse" << std::endl;

    TH1D *hRealTriggCounter = (TH1D*)fIn->Get("hRealTriggCounter");
    int nRealTrigg[nTriggBins];
    for (int it = 0; it < nTriggBins; it++) nRealTrigg[it] = hRealTriggCounter->GetBinContent(it+1);

    TH1D *hMassTrigg[nTriggBins];
    TH1D *hMassAssocPeak[nTriggBins][nAssocBins];
    TH1D *hMassAssocSide[nTriggBins][nAssocBins];

    TH2D *hCorrReal[nTriggBins][nAssocBins];
    TH2D *hCorrMassMass[nTriggBins][nAssocBins];
    TH2D *hCorrMassSide[nTriggBins][nAssocBins];
    TH2D *hCorrSideMass[nTriggBins][nAssocBins];
    TH2D *hCorrSideSide[nTriggBins][nAssocBins];
    for (int it = 0; it < nTriggBins; it++) {
        double tlow = triggPt[it];
        double tupp = triggPt[it+1];
        
        hMassTrigg[it] = (TH1D*)fIn->Get(Form("Masses/hPi0MassTrigg[%4.1f,%4.1f]",tlow,tupp));
        //hMassTrigg[it]->Scale(1./nEvent);
        //hMassTrigg[it]->Rebin();

        for (int ia = 0; ia < nAssocBins; ia++) {
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];

            if (tlow < aupp) continue;
            
            hMassAssocPeak[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocPeak[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            //hMassAssocPeak[it][ia]->Scale(1./nEvent);
            
            hMassAssocSide[it][ia] = (TH1D*)fIn->Get(Form("Masses/hPi0MassAssocSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp));
            //hMassAssocSide[it][ia]->Scale(1./nEvent);
            
            hCorrReal[it][ia] = (TH2D*)fIn->Get(Form("CorrFor/hCorrFor[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrReal[it][ia]->Rebin2D(4);
            hCorrMassMass[it][ia] = (TH2D*)fIn->Get(Form("CorrMassMass/hCorrMassMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrMassMass[it][ia]->Rebin2D(4);
            hCorrMassSide[it][ia] = (TH2D*)fIn->Get(Form("CorrMassSide/hCorrMassSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrMassSide[it][ia]->Rebin2D(4);
            hCorrSideMass[it][ia] = (TH2D*)fIn->Get(Form("CorrSideMass/hCorrSideMass[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrSideMass[it][ia]->Rebin2D(4);
            hCorrSideSide[it][ia] = (TH2D*)fIn->Get(Form("CorrSideSide/hCorrSideSide[%4.1f,%4.1f][%4.1f,%4.1f]",tlow,tupp,alow,aupp)); hCorrSideSide[it][ia]->Rebin2D(4);
        }
    }

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.58, 0.66, 0.78, 0.82);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.04); 

    TLegend *leg2 = new TLegend(0.58, 0.7, 0.78, 0.82);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.04); 

    // ----------------
    // |   Analysis   |
    // ----------------

    // 1. Fit mass distributions (to get parameters for normalisations)
    TF1 *fFitTrigg[nTriggBins];
    TF1 *fPeakTrigg[nTriggBins];
    TF1 *fBgTrigg[nTriggBins];

    TF1 *fPeakColored[nTriggBins];
    TF1 *fLeftSidebandColored[nTriggBins];
    TF1 *fRightSidebandColored[nTriggBins];
    
    double parTrigg[8];
    TCanvas *cMassTrigg[nTriggBins];
    for (int it = 0; it < nTriggBins; it++) {
        
        fFitTrigg[it] = new TF1("fFitTrigg", FitFunction, 20, 300, 8);
        fFitTrigg[it]->SetParameters(-1., 1., 1., 100., 100., 135., 10., 15.);
        fFitTrigg[it]->SetParLimits(5, 132, 138);
        fFitTrigg[it]->SetParLimits(6, 1., 15.);
        fFitTrigg[it]->SetParLimits(7, 2., 20.);
        fFitTrigg[it]->SetNpx(1000);
        
        fPeakTrigg[it] = new TF1("fPeakTrigg", FitPeak, 20, 300, 5);
        fPeakTrigg[it]->SetLineColor(kBlue);
        fPeakTrigg[it]->SetNpx(1000);
        
        fBgTrigg[it] = new TF1("fBgTrigg", FitBackground, 20, 300, 4);
        fBgTrigg[it]->SetLineColor(kBlack);
        fBgTrigg[it]->SetLineStyle(kDashed);
        fBgTrigg[it]->SetNpx(1000);
        
        hMassTrigg[it]->Fit("fFitTrigg", "Q0R+");
        fFitTrigg[it]->GetParameters(parTrigg);
        fBgTrigg[it]->SetParameters(parTrigg);
        fPeakTrigg[it]->SetParameters(&parTrigg[3]); 

        cMassTrigg[it] = new TCanvas(Form("cMassTrigg%d", it), Form("cMassTrigg%d", it), 600, 600);
        hMassTrigg[it]->SetTitle(Form("Trigger invariant mass, %.1f < p_{T,t} < %.1f GeV/c; M_{#gamma#gamma} (MeV/c^{2}); counts", triggPt[it], triggPt[it+1]));
        hMassTrigg[it]->GetXaxis()->SetRangeUser(0., 300.);
        hMassTrigg[it]->GetYaxis()->SetMaxDigits(3);
        
        fPeakColored[it] = (TF1*)fFitTrigg[it]->Clone("fPeakColored");
        fPeakColored[it]->SetFillColor(kBlue-10);
        fPeakColored[it]->GetXaxis()->SetRangeUser(110., 160.);
        fPeakColored[it]->SetFillStyle(1001);
        fLeftSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fLeftSidebandColored[it]->SetFillColor(kRed-10);
        fLeftSidebandColored[it]->GetXaxis()->SetRangeUser(40., 80.);
        fLeftSidebandColored[it]->SetFillStyle(1001);
        fRightSidebandColored[it] = (TF1*)fFitTrigg[it]->Clone("fLeftSidebandColored");
        fRightSidebandColored[it]->SetFillColor(kRed-10);
        fRightSidebandColored[it]->GetXaxis()->SetRangeUser(210., 280.);
        fRightSidebandColored[it]->SetFillStyle(1001);
        
        hMassTrigg[it]->Draw("HIST");
        fPeakColored[it]->Draw("SAME");
        fLeftSidebandColored[it]->Draw("SAME");
        fRightSidebandColored[it]->Draw("SAME");
        fFitTrigg[it]->Draw("SAME");
        fBgTrigg[it]->Draw("SAME");
        fPeakTrigg[it]->Draw("SAME");
        gPad->RedrawAxis();
    }

    TF1 *fFitAssoc[nTriggBins][nAssocBins];
    TF1 *fPeakAssoc[nTriggBins][nAssocBins];
    TF1 *fBgAssoc[nTriggBins][nAssocBins]; 
    
    TF1 *fPeakColoredAssoc[nTriggBins][nAssocBins];
    TF1 *fLeftSidebandColoredAssoc[nTriggBins][nAssocBins];
    TF1 *fRightSidebandColoredAssoc[nTriggBins][nAssocBins];

    double parAssoc[9];
    TCanvas *cMassAssocPeak[nTriggBins][nAssocBins];
    for (int it = 0; it < nTriggBins; it++) {
        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];
 
            if (tlow < aupp) continue;
            
            fFitAssoc[it][ia] = new TF1("fFitAssoc", FitFunction, 10, 300, 8);
            fFitAssoc[it][ia]->SetParameters(0., 0., 0., 0., 0., 135., 30., 50.);
            fFitAssoc[it][ia]->SetParLimits(5, 132, 138);
            fFitAssoc[it][ia]->SetParLimits(6, 2., 15.);
            fFitAssoc[it][ia]->SetParLimits(7, 4., 20.);

            fFitAssoc[it][ia]->SetNpx(1000);

            fPeakAssoc[it][ia] = new TF1("fPeakAssoc", FitPeak, 10, 300, 5);
            fPeakAssoc[it][ia]->SetLineColor(kBlue);
            fPeakAssoc[it][ia]->SetNpx(1000);

            fBgAssoc[it][ia]= new TF1("fBgAssoc", FitBackground, 10, 300, 3);
            fBgAssoc[it][ia]->SetLineColor(kBlack);
            fBgAssoc[it][ia]->SetLineStyle(kDashed);
            fBgAssoc[it][ia]->SetNpx(1000);

            hMassAssocPeak[it][ia]->Fit("fFitAssoc", "Q0R+");
            fFitAssoc[it][ia]->GetParameters(parAssoc);
            fBgAssoc[it][ia]->SetParameters(parAssoc);
            fPeakAssoc[it][ia]->SetParameters(&parAssoc[3]);
            
            cMassAssocPeak[it][ia] = new TCanvas(Form("cMassAssocPeak%d:%d",it,ia), Form("cMassAssocPeak%d:%d",it,ia), 600, 600);
            hMassAssocPeak[it][ia]->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; M_{#gamma#gamma}; counts", tlow, tupp, alow, aupp));
            hMassAssocPeak[it][ia]->GetXaxis()->SetRangeUser(0., 300.);
            hMassAssocPeak[it][ia]->GetYaxis()->SetMaxDigits(3);
            
            fPeakColoredAssoc[it][ia] = (TF1*)fFitAssoc[it][ia]->Clone("fPeakColoredAssoc");
            fPeakColoredAssoc[it][ia]->SetFillColor(kBlue-10);
            fPeakColoredAssoc[it][ia]->GetXaxis()->SetRangeUser(110., 160.);
            fPeakColoredAssoc[it][ia]->SetFillStyle(1001);
            fLeftSidebandColoredAssoc[it][ia] = (TF1*)fFitAssoc[it][ia]->Clone("fLeftSidebandColoredAssoc");
            fLeftSidebandColoredAssoc[it][ia]->SetFillColor(kRed-10);
            fLeftSidebandColoredAssoc[it][ia]->GetXaxis()->SetRangeUser(40., 80.);
            fLeftSidebandColoredAssoc[it][ia]->SetFillStyle(1001);
            fRightSidebandColoredAssoc[it][ia] = (TF1*)fFitAssoc[it][ia]->Clone("fLeftSidebandColoredAssoc");
            fRightSidebandColoredAssoc[it][ia]->SetFillColor(kRed-10);
            fRightSidebandColoredAssoc[it][ia]->GetXaxis()->SetRangeUser(210., 280.);
            fRightSidebandColoredAssoc[it][ia]->SetFillStyle(1001);
            
            hMassAssocPeak[it][ia]->Draw("HIST");
            fPeakColoredAssoc[it][ia]->Draw("SAME");
            fLeftSidebandColoredAssoc[it][ia]->Draw("SAME");
            fRightSidebandColoredAssoc[it][ia]->Draw("SAME");
            fFitAssoc[it][ia]->Draw("SAME");
            fBgAssoc[it][ia]->Draw("SAME");
            fPeakAssoc[it][ia]->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    // Integrate to get the normalisations
    double st[nTriggBins], bt[nTriggBins], sumt[nTriggBins], alphat[nTriggBins], betat[nTriggBins], bsidet[nTriggBins]; 
    double sa[nTriggBins][nAssocBins], ba[nTriggBins][nAssocBins], suma[nTriggBins][nAssocBins], alphaa[nTriggBins][nAssocBins], betaa[nTriggBins][nAssocBins], bsidea[nTriggBins][nAssocBins];
    double A[nTriggBins][nAssocBins], B[nTriggBins][nAssocBins];

    std::cout << "Number of triggers : " << std::endl;
    for (int it = 0; it < nTriggBins; it++) {
        sumt[it] = fFitTrigg[it]->Integral(110, 160);
        st[it] = fPeakTrigg[it]->Integral(110, 160);
        bt[it] = fBgTrigg[it]->Integral(110, 160);
        bsidet[it] = fBgTrigg[it]->Integral(40, 80) + fBgTrigg[it]->Integral(210, 280);
        alphat[it] = st[it]/sumt[it];
        betat[it] = bt[it]/sumt[it];
        
        for (int ia = 0; ia < nAssocBins; ia++) {
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];
 
            if (tlow < aupp) continue;
            
            suma[it][ia] = fFitAssoc[it][ia]->Integral(110, 160);
            sa[it][ia] = fPeakAssoc[it][ia]->Integral(110, 160);
            ba[it][ia] = fBgAssoc[it][ia]->Integral(110, 160);
            bsidea[it][ia] = fBgAssoc[it][ia]->Integral(40, 80) + fBgAssoc[it][ia]->Integral(210, 280);
            alphaa[it][ia] = sa[it][ia]/suma[it][ia];
            betaa[it][ia] = ba[it][ia]/suma[it][ia];
            A[it][ia] = bt[it]/bsidet[it];
            B[it][ia] = ba[it][ia]/bsidea[it][ia];
        }
        
        std::cout << "\tbin [ " << triggPt[it] << " " << triggPt[it+1] << " ] : \treal=" << nRealTrigg[it] << "\treconst=" << st[it] << std::endl;
        
    }

    TCanvas *cCorrFuncs[nTriggBins][nAssocBins]; 
    TCanvas *cCorr[nTriggBins][nAssocBins]; 
    TCanvas *cRatio[nTriggBins][nAssocBins];
    TH1D *hCorrRealProj[nTriggBins][nAssocBins];
    TH1D *hCorrMassMassProj[nTriggBins][nAssocBins];
    TH1D *hCorrMassSideProj[nTriggBins][nAssocBins];
    TH1D *hCorrSideMassProj[nTriggBins][nAssocBins];
    TH1D *hCorrSideSideProj[nTriggBins][nAssocBins];
    TH1D *hCorr[nTriggBins][nAssocBins];
    TH1D *hRatio[nTriggBins][nAssocBins];
    
    // 2. Normalise correlation functions   
    for (int it = 0; it < nTriggBins; it++) {
        for (int ia = 0; ia < nAssocBins; ia++) {        
            double tlow = triggPt[it];
            double tupp = triggPt[it+1];
            double alow = assocPt[ia];
            double aupp = assocPt[ia+1];
 
            if (tlow < aupp) continue;
            
            double nPairReal = hCorrReal[it][ia]->GetEntries();
            double nPairMassMass = hCorrMassMass[it][ia]->GetEntries();
            double nPairMassSide = hCorrMassSide[it][ia]->GetEntries();
            double nPairSideMass = hCorrSideMass[it][ia]->GetEntries();
            double nPairSideSide = hCorrSideSide[it][ia]->GetEntries();
            
            hCorrRealProj[it][ia] = hCorrReal[it][ia]->ProjectionX();
            //hCorrRealProj[it][ia]->Scale(1.0/nPairMassMass);
            hCorrRealProj[it][ia]->SetLineColor(kRed);

            hCorrMassMassProj[it][ia] = hCorrMassMass[it][ia]->ProjectionX();
            //hCorrMassMassProj[it][ia]->Scale(1.0/nPairMassMass);
            //hCorrMassMassProj[it][ia]->Scale(1/(alphat[it]*alphaa[it][ia]));
            hCorrMassMassProj[it][ia]->SetLineColor(kBlue);
            hCorrMassMassProj[it][ia]->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; dN/d#Delta#phi", triggPt[it], triggPt[it+1], assocPt[ia], assocPt[ia+1]));
            
            hCorrMassSideProj[it][ia] = hCorrMassSide[it][ia]->ProjectionX();
            //hCorrMassSideProj[it][ia]->Scale(1.0/nPairMassSide);
            //hCorrMassSideProj[it][ia]->Scale(betaa[it][ia]/(alphat[it]*alphaa[it][ia]));
            hCorrMassSideProj[it][ia]->Scale(B[it][ia]);
            hCorrMassSideProj[it][ia]->SetLineColor(kRed);
            
            hCorrSideMassProj[it][ia] = hCorrSideMass[it][ia]->ProjectionX();
            //hCorrSideMassProj[it][ia]->Scale(1.0/nPairSideMass);
            //hCorrSideMassProj[it][ia]->Scale(betat[it]/(alphat[it]*alphaa[it][ia]));
            hCorrSideMassProj[it][ia]->Scale(A[it][ia]);
            hCorrSideMassProj[it][ia]->SetLineColor(kBlack);
            
            hCorrSideSideProj[it][ia] = hCorrSideSide[it][ia]->ProjectionX();
            //hCorrSideSideProj[it][ia]->Scale(1.0/nPairSideSide);
            //hCorrSideSideProj[it][ia]->Scale((betat[it]*betaa[it][ia])/(alphat[it]*alphaa[it][ia]));
            hCorrSideSideProj[it][ia]->Scale((A[it][ia]*B[it][ia]));
            hCorrSideSideProj[it][ia]->SetLineColor(kGray);
            
            if (it==0 && ia==0) {
                leg1->AddEntry(hCorrMassMassProj[it][ia], "f_{mass,mass}", "l");
                leg1->AddEntry(hCorrMassSideProj[it][ia], "f_{mass,side}", "l");
                leg1->AddEntry(hCorrSideMassProj[it][ia], "f_{side,mass}", "l");
                leg1->AddEntry(hCorrSideSideProj[it][ia], "f_{side,side}", "l");
            }

            /**cCorrFuncs[it][ia] = new TCanvas(Form("cCorrFuncs%d:%d",it,ia), Form("cCorrFuncs%d:%d",it,ia), 600, 600);
            hCorrMassMassProj[it][ia]->Draw("HIST");
            hCorrMassSideProj[it][ia]->Draw("HIST SAME");
            hCorrSideMassProj[it][ia]->Draw("HIST SAME");
            hCorrSideSideProj[it][ia]->Draw("HIST SAME");
            leg1->Draw("SAME");**/
    
            // 3. Get the final correlation function with sideband correction
            cCorr[it][ia] = new TCanvas(Form("cCorr%d:%d",it,ia), Form("cCorr%d:%d",it,ia), 600, 600);

            hCorr[it][ia] = (TH1D*)hCorrMassMassProj[it][ia]->Clone(Form("hCorr%d:%d",it,ia));
            hCorr[it][ia]->SetLineColor(kBlack);
            hCorr[it][ia]->Add(hCorrMassSideProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideMassProj[it][ia], -1);
            hCorr[it][ia]->Add(hCorrSideSideProj[it][ia]);
            hCorr[it][ia]->Scale(1./(brGammaCh*brGammaCh)); // Correct for branching ratio
            hCorr[it][ia]->Scale(1./(0.95*0.95)); // Correct for missing photon pair at limit of acceptance

            if (it==0 && ia==0) {
                leg2->AddEntry(hCorrMassMassProj[it][ia], "f_{mass,mass}", "l");
                leg2->AddEntry(hCorr[it][ia], "f_{SS}", "l");
                leg2->AddEntry(hCorrRealProj[it][ia], "f_{real}", "l");
            }

            hCorrMassMassProj[it][ia]->Draw("HIST");
            hCorr[it][ia]->Draw("HIST SAME");
            hCorrRealProj[it][ia]->Draw("HIST SAME");
            leg2->Draw("SAME");

            cRatio[it][ia] = new TCanvas(Form("cRatio%d:%d",it,ia), Form("cRatio%d:%d",it,ia), 600, 600);
            
            hRatio[it][ia] = (TH1D*)hCorr[it][ia]->Clone(Form("hRatio%d:%d",it,ia));
            hRatio[it][ia]->Divide(hCorrRealProj[it][ia]);
            hRatio[it][ia]->GetYaxis()->SetRangeUser(0., 2.);
            hRatio[it][ia]->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; f_{SS}/f_{real}", triggPt[it], triggPt[it+1], assocPt[ia], assocPt[ia+1]));
            hRatio[it][ia]->Fit("pol1");
            hRatio[it][ia]->Draw();
        }
    }
}


// -----------------
// |   Functions   |
// -----------------

double FitPeak(double *x, double *p)
{   
    double c1 = p[0];
    double c2 = p[1];
    double mu = p[2];
    double sigma1 = p[3];
    double sigma2 = p[4];

    return c1*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma1*sigma1)) + c2*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma2*sigma2));
}

double FitBackground(double *x, double *p)
{
    double a = p[0];
    double b = p[1];
    double c = p[2];

    return a*x[0]*x[0] + b*x[0] + c;
}

double FitFunction(double *x, double *p)
{
    return FitBackground(x, p) + FitPeak(x, &p[3]);
}
