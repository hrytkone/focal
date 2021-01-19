const int nTriggBins = 1;
double  triggPt[nTriggBins+1] = {5.0, 15.0};

const int nAssocBins = 1;
double  assocPt[nAssocBins+1] = {2.0, 5.0};

double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);

void AnalyseCorrelations2(TString sInputName = "output.root")
{

    gStyle->SetOptStat(0);

    TFile *fIn = TFile::Open(sInputName);

    // ------------------
    // |   Histograms   |
    // ------------------
    TH1D *hCounter = (TH1D*)fIn->Get("hCounter");
    int nEvent = hCounter->GetBinContent(1);
    std::cout << "Input file contains " << nEvent << " events, proceed to analyse" << std::endl;

    TH2D *hCorrReal[nTriggBins][nAssocBins];

    TH2D *hCorrMassMass[nTriggBins][nAssocBins];
    TH2D *hCorrMassSide[nTriggBins][nAssocBins];
    TH2D *hCorrSideMass[nTriggBins][nAssocBins];
    TH2D *hCorrSideSide[nTriggBins][nAssocBins];
    for (int it = 0; it < nTriggBins; it++) {
        for (int ia = 0; ia < nAssocBins; ia++) {
            hCorrReal[it][ia] = (TH2D*)fIn->Get(Form("hCorrFor%d:%d", it, ia)); hCorrReal[it][ia]->Rebin2D();
            hCorrMassMass[it][ia] = (TH2D*)fIn->Get(Form("hCorrMassMass%d:%d", it, ia)); hCorrMassMass[it][ia]->Rebin2D();
            hCorrMassSide[it][ia] = (TH2D*)fIn->Get(Form("hCorrMassSide%d:%d", it, ia)); hCorrMassSide[it][ia]->Rebin2D();
            hCorrSideMass[it][ia] = (TH2D*)fIn->Get(Form("hCorrSideMass%d:%d", it, ia)); hCorrSideMass[it][ia]->Rebin2D();
            hCorrSideSide[it][ia] = (TH2D*)fIn->Get(Form("hCorrSideSide%d:%d", it, ia)); hCorrSideSide[it][ia]->Rebin2D();
        }
    }

    double nTriggEvent = hCounter->GetBinContent(4);
    TH1D *hMassMassTrigg[nTriggBins];
    TH1D *hMassSideTrigg[nTriggBins];
    TH1D *hSideMassTrigg[nTriggBins];
    TH1D *hSideSideTrigg[nTriggBins];
    TH1D *hMassMassAssoc[nAssocBins];
    TH1D *hMassSideAssoc[nAssocBins];
    TH1D *hSideMassAssoc[nAssocBins];
    TH1D *hSideSideAssoc[nAssocBins];
    for (int it = 0; it < nTriggBins; it++) {
        hMassMassTrigg[it] = (TH1D*)fIn->Get(Form("hPi0MassMassTrigg%d", it)); //hMassMassTrigg[it]->Rebin();
        hMassMassTrigg[it]->Scale(1./nTriggEvent);
        hMassSideTrigg[it] = (TH1D*)fIn->Get(Form("hPi0MassSideTrigg%d", it)); //hMassSideTrigg[it]->Rebin();
        hMassSideTrigg[it]->Scale(1./nTriggEvent);
        hSideMassTrigg[it] = (TH1D*)fIn->Get(Form("hPi0SideMassTrigg%d", it)); //hSideMassTrigg[it]->Rebin();
        hSideMassTrigg[it]->Scale(1./nTriggEvent);
        hSideSideTrigg[it] = (TH1D*)fIn->Get(Form("hPi0SideSideTrigg%d", it)); //hSideSideTrigg[it]->Rebin();
        hSideSideTrigg[it]->Scale(1./nTriggEvent);
    }

    for (int ia = 0; ia < nAssocBins; ia++) {
        hMassMassAssoc[ia] = (TH1D*)fIn->Get(Form("hPi0MassMassAssoc%d", ia)); //hMassMassAssoc[ia]->Rebin();
        hMassMassAssoc[ia]->Scale(1./nTriggEvent);
        hMassSideAssoc[ia] = (TH1D*)fIn->Get(Form("hPi0MassSideAssoc%d", ia)); //hMassSideAssoc[ia]->Rebin();
        hMassSideAssoc[ia]->Scale(1./nTriggEvent);
        hSideMassAssoc[ia] = (TH1D*)fIn->Get(Form("hPi0SideMassAssoc%d", ia)); //hSideMassAssoc[ia]->Rebin();
        hSideMassAssoc[ia]->Scale(1./nTriggEvent);
        hSideSideAssoc[ia] = (TH1D*)fIn->Get(Form("hPi0SideSideAssoc%d", ia)); //hSideSideAssoc[ia]->Rebin();
        hSideSideAssoc[ia]->Scale(1./nTriggEvent);
    }

    // ---------------
    // |   Legends   |
    // ---------------
    TLegend *leg1 = new TLegend(0.58, 0.6, 0.78, 0.85);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.05); 

    TLegend *leg2 = new TLegend(0.58, 0.66, 0.78, 0.85);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.05); 

    // ----------------
    // |   Analysis   |
    // ----------------
    // Data from Pythia pi0s
    std::cout << "\nReal pi0 statistics: " << std::endl;
    std::cout << "\tNumber of events with triggers: " << hCounter->GetBinContent(4) << std::endl;
    std::cout << "\tNumber of triggers: " << hCounter->GetBinContent(5) << std::endl;
    std::cout << "\tTriggers per event (all): " << hCounter->GetBinContent(5)/nEvent << std::endl;
    std::cout << "\tTriggers per event (with trigger): " << hCounter->GetBinContent(5)/hCounter->GetBinContent(4) << std::endl;
    std::cout << "\tAssociated per trigger: " << hCounter->GetBinContent(6)/hCounter->GetBinContent(5) << "\n" << std::endl;
    

    // 1. Fit mass distributions (to get parameters for normalisations)
    TF1 *fFitTrigg = new TF1("fFitTrigg", FitFunction, 20, 300, 6);
    //fFitTrigg->SetParameters(0., 0., 0., 0., 0., 135., 30., 60.);
    //fFitTrigg->SetParameters(0., 0., 0.1, 1.0, 135.0);
    fFitTrigg->SetParLimits(4, 5., 20.);
    fFitTrigg->SetParLimits(5, 130, 140);
    fFitTrigg->SetNpx(1000);
    
    TF1 *fPeakTrigg = new TF1("fPeakTrigg", FitPeak, 20, 300, 6);
    fPeakTrigg->SetLineColor(kBlue);
    fPeakTrigg->SetNpx(1000);

    TF1 *fBgTrigg = new TF1("fBgTrigg", FitBackground, 20, 300, 3);
    fBgTrigg->SetLineColor(kBlack);
    fBgTrigg->SetLineStyle(kDashed);
    fBgTrigg->SetNpx(1000);

    double parTrigg[9];
    TCanvas *cMassTrigg = new TCanvas("cMassTrigg", "cMassTrigg");
    for (int it = 0; it < nTriggBins; it++) {
        hMassMassTrigg[it]->Fit("fFitTrigg", "R");
        fFitTrigg->GetParameters(parTrigg);
        fBgTrigg->SetParameters(parTrigg);
        fPeakTrigg->SetParameters(&parTrigg[3]);

        hMassMassTrigg[it]->SetTitle(Form("Trigger invariant mass, %.1f < p_{T,t} < %.1f GeV/c", triggPt[it], triggPt[it+1]));

        hMassMassTrigg[it]->Draw("HIST");
        fFitTrigg->Draw("SAME");
        fBgTrigg->Draw("SAME");
        fPeakTrigg->Draw("SAME");
    }

    TF1 *fFitAssoc = new TF1("fFitAssoc", FitFunction, 10, 300, 6);
    //fFitAssoc->SetParameters(-0.1, 0.1, 0.1, 1., 1., 135., 10., 15.);
    fFitAssoc->SetParLimits(4, 5., 30.);
    fFitAssoc->SetParLimits(5, 130, 140);
    fFitAssoc->SetNpx(1000);

    TF1 *fPeakAssoc = new TF1("fPeakAssoc", FitPeak, 10, 300, 6);
    fPeakAssoc->SetLineColor(kBlue);
    fPeakAssoc->SetNpx(1000);

    TF1 *fBgAssoc = new TF1("fBgAssoc", FitBackground, 10, 300, 3);
    fBgAssoc->SetLineColor(kBlack);
    fBgAssoc->SetLineStyle(kDashed);
    fBgAssoc->SetNpx(1000);

    double parAssoc[9];
    TCanvas *cMassAssoc = new TCanvas("cMassAssoc", "cMassAssoc");
    for (int ia = 0; ia < nAssocBins; ia++) {
        hMassMassAssoc[ia]->Fit("fFitAssoc", "R");
        fFitAssoc->GetParameters(parAssoc);
        fBgAssoc->SetParameters(parAssoc);
        fPeakAssoc->SetParameters(&parAssoc[3]);
        
        hMassMassAssoc[ia]->SetTitle(Form("Associated invariant mass, %.1f < p_{T,a} < %.1f GeV/c", assocPt[ia], assocPt[ia+1]));
        
        hMassMassAssoc[ia]->Draw("HIST");
        fFitAssoc->Draw("SAME");
        fBgAssoc->Draw("SAME");
        fPeakAssoc->Draw("SAME");
    }
    
    //TF1 *fFitTriggMassSide = new TF1("fFitTriggMassSide", FitFunction, 20, 300, 8);
    TF1 *fFitTriggMassSide = new TF1("fFitTriggMassSide", FitFunction, 20, 300, 6);
    //fFitTriggMassSide->SetParameters(0., 0., 0., 0., 0., 135., 1., 1.);
    fFitTriggMassSide->SetParLimits(4, 5., 20.);
    fFitTriggMassSide->SetParLimits(5, 130, 140);
    fFitTriggMassSide->SetNpx(1000);
    
    TF1 *fPeakTriggMassSide = new TF1("fPeakTriggMassSide", FitPeak, 20, 300, 6);
    fPeakTriggMassSide->SetLineColor(kBlue);
    fPeakTriggMassSide->SetNpx(1000);

    TF1 *fBgTriggMassSide = new TF1("fBgTriggMassSide", FitBackground, 20, 300, 3);
    fBgTriggMassSide->SetLineColor(kBlack);
    fBgTriggMassSide->SetLineStyle(kDashed);
    fBgTriggMassSide->SetNpx(1000);

    double parTriggMassSide[9];
    TCanvas *cMassTriggMassSide = new TCanvas("cMassTriggMassSide", "cMassTriggMassSide");
    for (int it = 0; it < nTriggBins; it++) {
        hMassSideTrigg[it]->Fit("fFitTriggMassSide", "R");
        fFitTriggMassSide->GetParameters(parTriggMassSide);
        fBgTriggMassSide->SetParameters(parTriggMassSide);
        fPeakTriggMassSide->SetParameters(&parTriggMassSide[3]);

        hMassSideTrigg[it]->SetTitle(Form("Trigger invariant mass (mass-side), %.1f < p_{T,t} < %.1f GeV/c", triggPt[it], triggPt[it+1]));

        hMassSideTrigg[it]->Draw("HIST");
        fFitTriggMassSide->Draw("SAME");
        fBgTriggMassSide->Draw("SAME");
        fPeakTriggMassSide->Draw("SAME");
    }
   
    //TF1 *fFitAssocSideMass = new TF1("fFitAssocSideMass", FitFunction, 10, 300, 8);
    TF1 *fFitAssocSideMass = new TF1("fFitAssocSideMass", FitFunction, 10, 300, 6);
    fFitAssocSideMass->SetParLimits(4, 5, 20);
    fFitAssocSideMass->SetParLimits(5, 130, 140);
    //fFitAssocSideMass->SetParameters(-0.1, 0.1, 0., 0., 0., 135., 30., 60.);
    fFitAssocSideMass->SetNpx(1000);

    TF1 *fPeakAssocSideMass = new TF1("fPeakAssocSideMass", FitPeak, 10, 300, 6);
    fPeakAssocSideMass->SetLineColor(kBlue);
    fPeakAssocSideMass->SetNpx(1000);

    TF1 *fBgAssocSideMass = new TF1("fBgAssocSideMass", FitBackground, 10, 300, 3);
    fBgAssocSideMass->SetLineColor(kBlack);
    fBgAssocSideMass->SetLineStyle(kDashed);
    fBgAssocSideMass->SetNpx(1000);

    double parAssocSideMass[9];
    TCanvas *cMassAssocSideMass = new TCanvas("cMassAssocSideMass", "cMassAssocSideMass");
    for (int ia = 0; ia < nAssocBins; ia++) {
        hSideMassAssoc[ia]->Fit("fFitAssocSideMass", "R");
        fFitAssocSideMass->GetParameters(parAssocSideMass);
        fBgAssocSideMass->SetParameters(parAssocSideMass);
        fPeakAssocSideMass->SetParameters(&parAssocSideMass[3]);
        
        hSideMassAssoc[ia]->SetTitle(Form("Associated invariant mass (side-mass), %.1f < p_{T,a} < %.1f GeV/c", assocPt[ia], assocPt[ia+1]));
        
        hSideMassAssoc[ia]->Draw("HIST");
        fFitAssocSideMass->Draw("SAME");
        fBgAssocSideMass->Draw("SAME");
        fPeakAssocSideMass->Draw("SAME");
    }

    // Integrate to get the normalisation
    double a1 = fPeakTrigg->Integral(120, 150);
    double b1 = fBgTrigg->Integral(120, 150);
    double c1 = fPeakAssoc->Integral(120, 150);
    double d1 = fBgAssoc->Integral(120, 150);

    double a2 = fPeakTriggMassSide->Integral(120, 150);
    double b2 = fBgTriggMassSide->Integral(120, 150);
    double f2 = hMassSideAssoc[0]->Integral(40, 80) + hMassSideAssoc[0]->Integral(210, 280);

    double e3 = hSideMassTrigg[0]->Integral(40, 80) + hSideMassTrigg[0]->Integral(210, 280);
    double c3 = fPeakAssocSideMass->Integral(120, 150);
    double d3 = fBgAssocSideMass->Integral(120, 150);
   
    double e4 = hSideSideTrigg[0]->Integral(40, 80) + hSideSideTrigg[0]->Integral(210, 280);
    double f4 = hSideSideAssoc[0]->Integral(40, 80) + hSideSideAssoc[0]->Integral(210, 280);

    double n1 = (a1*d1)/(a2*f2);
    double n2 = (b1*c1)/(e3*c3);
    double n3 = (b1*d1)/(e4*f4)*((c1*d3)/(d1*c3) + (a1*b2)*(a2*b1) - 1);

    std::cout << "\nMass-mass : " << std::endl;
    std::cout << "\ta1 : " << a1 << std::endl;
    std::cout << "\tb1 : " << b1 << std::endl;
    std::cout << "\tc1 : " << c1 << std::endl;
    std::cout << "\td1 : " << d1 << std::endl;
    std::cout << "\nMass-side : " << std::endl;
    std::cout << "\ta2 : " << a2 << std::endl;
    std::cout << "\tb2 : " << b2 << std::endl;
    std::cout << "\tf2 : " << f2 << std::endl;
    std::cout << "\nSide-mass : " << std::endl;
    std::cout << "\te3 : " << e3 << std::endl;
    std::cout << "\tc3 : " << c3 << std::endl;
    std::cout << "\td3 : " << d3 << std::endl;
    std::cout << "\nSide-side : " << std::endl;
    std::cout << "\te4 : " << e4 << std::endl;
    std::cout << "\tf4 : " << f4 << std::endl;
    std::cout << "\nNormalisations : " << std::endl;
    std::cout << "\tn1 : " << n1 << std::endl;
    std::cout << "\tn2 : " << n2 << std::endl;
    std::cout << "\tn3 : " << n3 << std::endl;
    
    // 2. Normalise correlation functions   
    for (int it = 0; it < nTriggBins; it++) {
        for (int ia = 0; ia < nAssocBins; ia++) {        

            TCanvas *cCorrFuncs = new TCanvas("cCorrFuncs", "cCorrFuncs");
            
            double nPairReal = hCorrReal[it][ia]->GetEntries();
            double nPairMassMass = hCorrMassMass[it][ia]->GetEntries();
            double nPairMassSide = hCorrMassSide[it][ia]->GetEntries();
            double nPairSideMass = hCorrSideMass[it][ia]->GetEntries();
            double nPairSideSide = hCorrSideSide[it][ia]->GetEntries();
            
            TH1D *hCorrRealProj = hCorrReal[it][ia]->ProjectionX();
            //hCorrRealProj->Scale(1.0/nPairMassMass);
            hCorrRealProj->SetLineColor(kRed);

            TH1D *hCorrMassMassProj = hCorrMassMass[it][ia]->ProjectionX();
            //hCorrMassMassProj->Scale(1.0/nPairMassMass);
            hCorrMassMassProj->SetLineColor(kBlue);
            hCorrMassMassProj->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; 1/N_{pair}dN/d#Delta#phi", triggPt[it], triggPt[it+1], assocPt[ia], assocPt[ia+1]));
            
            TH1D *hCorrMassSideProj = hCorrMassSide[it][ia]->ProjectionX();
            //hCorrMassSideProj->Scale(1.0/nPairMassSide);
            hCorrMassSideProj->Scale(n1);
            hCorrMassSideProj->SetLineColor(kRed);
            
            TH1D *hCorrSideMassProj = hCorrSideMass[it][ia]->ProjectionX();
            //hCorrSideMassProj->Scale(1.0/nPairSideMass);
            hCorrSideMassProj->Scale(n2);
            hCorrSideMassProj->SetLineColor(kBlack);
            
            TH1D *hCorrSideSideProj = hCorrSideSide[it][ia]->ProjectionX();
            //hCorrSideSideProj->Scale(1.0/nPairSideSide);
            hCorrSideSideProj->Scale(n3);
            hCorrSideSideProj->SetLineColor(kGray);

            leg1->AddEntry(hCorrMassMassProj, "f_{mass,mass}", "l");
            leg1->AddEntry(hCorrMassSideProj, "f_{mass,side}", "l");
            leg1->AddEntry(hCorrSideMassProj, "f_{side,mass}", "l");
            leg1->AddEntry(hCorrSideSideProj, "f_{side,side}", "l");

            hCorrMassMassProj->Draw("HIST");
            hCorrMassSideProj->Draw("HIST SAME");
            hCorrSideMassProj->Draw("HIST SAME");
            hCorrSideSideProj->Draw("HIST SAME");
            leg1->Draw("SAME");
    
            // 3. Get the final correlation function with sideband correction
            TCanvas *cCorr = new TCanvas("cCorr", "cCorr");
            TH1D *hCorr = (TH1D*)hCorrMassMassProj->Clone("hCorr");
            hCorr->SetLineColor(kBlack);
            hCorr->Add(hCorrMassSideProj, -1);
            hCorr->Add(hCorrSideMassProj, -1);
            hCorr->Add(hCorrSideSideProj);

            leg2->AddEntry(hCorrMassMassProj, "f_{mass,mass}", "l");
            leg2->AddEntry(hCorr, "f_{signal,signal}", "l");
            leg2->AddEntry(hCorrRealProj, "f_{real}", "l");

            hCorrMassMassProj->Draw("HIST");
            hCorr->Draw("HIST SAME");
            hCorrRealProj->Draw("HIST SAME");
            leg2->Draw("SAME");

            TCanvas *cRatio = new TCanvas("cRatio", "cRatio");
            TH1D *hRatio = (TH1D*)hCorr->Clone("hRatio");
            hRatio->Divide(hCorrRealProj);
            hRatio->SetTitle("f_{mass,mass}/f_{real}; #Delta#phi;");
            hRatio->Draw();

            TCanvas *cSideMassVsMassSide = new TCanvas("cSideMassVsMassSide");
            TH1D *hSideMassVsMassSide = (TH1D*)hCorrSideMassProj->Clone("hSideMassVsMassSide");
            hSideMassVsMassSide->Divide(hCorrMassSideProj);
            hSideMassVsMassSide->SetTitle("f_{side,mass}/f_{mass,side}; #Delta#phi;");
            hSideMassVsMassSide->Draw();

            TCanvas *cAllVsMassMass = new TCanvas("cAllVsMassMass");
            TH1D *hAllVsMassMass = (TH1D*)hCorr->Clone("hAllVsMassMass");
            hAllVsMassMass->Divide(hCorrMassMassProj);
            hAllVsMassMass->SetTitle("f_{signal,signal}/f_{mass,mass}; #Delta#phi;");
            hAllVsMassMass->Draw();

        }
    }

}


// -----------------
// |   Functions   |
// -----------------

/**double FitPeak(double *x, double *p)
{   
    double c1 = p[0];
    double c2 = p[1];
    double mu = p[2];
    double sigma1 = p[3];
    double sigma2 = p[4];

    return c1*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma1*sigma1)) + c2*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2*sigma2*sigma2));
}**/

double FitPeak(double *x, double *p) 
{
    double b = p[0];
    double c = p[1];
    double m = p[2];
    return (0.5*b*c/TMath::Pi()) / TMath::Max(1.e-10, (x[0]-m)*(x[0]-m)+ .25*c*c);
}

/**double FitPeak(double *x, double *p)
{
    double b = p[0];
    double mu = p[1];

    return TMath::Exp(-TMath::Abs((x[0]-mu)/b))/(2*b);
}**/

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

