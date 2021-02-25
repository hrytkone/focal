const int nTriggBins = 1;
double  triggPt[nTriggBins+1] = {5.0, 15.0};

const int nAssocBins = 1;
double  assocPt[nAssocBins+1] = {2.0, 5.0};

double FitPeak(double *x, double *p);
double FitBackground(double *x, double *p);
double FitFunction(double *x, double *p);

void AnalyseCorrelations(TString sInputName = "output.root")
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
            hCorrReal[it][ia] = (TH2D*)fIn->Get(Form("hCorrFor%d:%d", it, ia)); hCorrReal[it][ia]->Rebin();
            hCorrMassMass[it][ia] = (TH2D*)fIn->Get(Form("hCorrMassMass%d:%d", it, ia)); hCorrMassMass[it][ia]->Rebin();
            hCorrMassSide[it][ia] = (TH2D*)fIn->Get(Form("hCorrMassSide%d:%d", it, ia)); hCorrMassSide[it][ia]->Rebin();
            hCorrSideMass[it][ia] = (TH2D*)fIn->Get(Form("hCorrSideMass%d:%d", it, ia)); hCorrSideMass[it][ia]->Rebin();
            hCorrSideSide[it][ia] = (TH2D*)fIn->Get(Form("hCorrSideSide%d:%d", it, ia)); hCorrSideSide[it][ia]->Rebin();
        }
    }

    TH1D *hMassTrigg[nTriggBins];
    TH1D *hMassAssoc[nAssocBins];
    for (int it = 0; it < nTriggBins; it++) {
        hMassTrigg[it] = (TH1D*)fIn->Get(Form("hPi0MassTrigg%d", it));
        hMassTrigg[it]->Scale(1./nEvent);
        hMassTrigg[it]->Rebin();
    }

    for (int ia = 0; ia < nAssocBins; ia++) {
        hMassAssoc[ia] = (TH1D*)fIn->Get(Form("hPi0MassAssoc%d", ia));
        hMassAssoc[ia]->Scale(1./nEvent);
        hMassAssoc[ia]->Rebin();
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
    TF1 *fFitTrigg = new TF1("fFitTrigg", FitFunction, 20, 300, 9);
    fFitTrigg->SetParameters(0., 0., 0., 0., 0., 135., 135., 30., 60.);
    //fFitTrigg->SetParLimits(1, 0., 1000.);
    fFitTrigg->SetParLimits(6, 130., 140.);
    fFitTrigg->SetParLimits(8, -10., 10.);
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
        hMassTrigg[it]->Fit("fFitTrigg", "R");
        fFitTrigg->GetParameters(parTrigg);
        fBgTrigg->SetParameters(parTrigg);
        fPeakTrigg->SetParameters(&parTrigg[3]);

        hMassTrigg[it]->SetTitle(Form("Trigger invariant mass, %.1f < p_{T,t} < %.1f GeV/c", triggPt[it], triggPt[it+1]));

        hMassTrigg[it]->Draw("HIST");
        fFitTrigg->Draw("SAME");
        fBgTrigg->Draw("SAME");
        fPeakTrigg->Draw("SAME");
    }

    TF1 *fFitAssoc = new TF1("fFitAssoc", FitFunction, 10, 300, 9);
    fFitAssoc->SetParameters(0., 0., 0., 0., 0., 135., 135., 30., 60.);
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
        hMassAssoc[ia]->Fit("fFitAssoc", "R");
        fFitAssoc->GetParameters(parAssoc);
        fBgAssoc->SetParameters(parAssoc);
        fPeakAssoc->SetParameters(&parAssoc[3]);
        
        hMassAssoc[ia]->SetTitle(Form("Associated invariant mass, %.1f < p_{T,a} < %.1f GeV/c", assocPt[ia], assocPt[ia+1]));
        
        hMassAssoc[ia]->Draw("HIST");
        fFitAssoc->Draw("SAME");
        fBgAssoc->Draw("SAME");
        fPeakAssoc->Draw("SAME");
    }

    // Integrate to get the normalisations
    double st, bt, sumt, sa, ba, suma;

    sumt = fFitTrigg->Integral(110, 170);
    st = fPeakTrigg->Integral(110, 170);
    bt = fBgTrigg->Integral(110, 170);
    double btside = fBgTrigg->Integral(55, 75) + fBgTrigg->Integral(210, 280);

    suma = fFitAssoc->Integral(110, 170);
    sa = fPeakAssoc->Integral(110, 170);
    ba = fBgAssoc->Integral(110, 170);
    double baside = fBgAssoc->Integral(55, 75) + fBgAssoc->Integral(210, 280);

    double alphat, alphaa, betat, betaa;
    alphat = st/sumt;
    alphaa = sa/suma;
    betat = bt/sumt;
    betaa = ba/suma;

    std::cout << "\nTrigger : " << std::endl;
    std::cout << "\tSignal : " << st << std::endl;
    std::cout << "\tBg (peak) : " << bt << std::endl;
    std::cout << "\tBg (side) : " << btside << std::endl;
    std::cout << "\tTotal : " << sumt << std::endl;
    std::cout << "\n\talpha : " << alphat << std::endl;
    std::cout << "\tbeta : " << betat << std::endl;
    std::cout << "\nAssoc : " << std::endl;
    std::cout << "\tSignal : " << sa << std::endl;
    std::cout << "\tBg (peak) : " << ba << std::endl;
    std::cout << "\tBg (side) : " << baside << std::endl;
    std::cout << "\tTotal : " << suma << std::endl;
    std::cout << "\n\talpha : " << alphaa << std::endl;
    std::cout << "\tbeta : " << betaa << std::endl;
 

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
            hCorrRealProj->Scale(1.0/nPairMassMass);
            hCorrRealProj->SetLineColor(kRed);

            TH1D *hCorrMassMassProj = hCorrMassMass[it][ia]->ProjectionX();
            hCorrMassMassProj->Scale(1.0/nPairMassMass);
            hCorrMassMassProj->Scale(1/(alphat*alphaa));
            //hCorrMassMassProj->Scale(1/st*sa);
            hCorrMassMassProj->SetLineColor(kBlue);
            hCorrMassMassProj->SetTitle(Form("%.1f < p_{T,t} < %.1f GeV/c, %.1f < p_{T,a} < %.1f GeV/c; #Delta#phi; 1/N_{pair}dN/d#Delta#phi", triggPt[it], triggPt[it+1], assocPt[ia], assocPt[ia+1]));
            
            TH1D *hCorrMassSideProj = hCorrMassSide[it][ia]->ProjectionX();
            hCorrMassSideProj->Scale(1.0/nPairMassSide);
            hCorrMassSideProj->Scale(betaa/(alphat*alphaa));
            //hCorrMassSideProj->Scale(ba/baside);
            hCorrMassSideProj->SetLineColor(kRed);
            
            TH1D *hCorrSideMassProj = hCorrSideMass[it][ia]->ProjectionX();
            hCorrSideMassProj->Scale(1.0/nPairSideMass);
            hCorrSideMassProj->Scale(betat/(alphat*alphaa));
            //hCorrSideMassProj->Scale(bt/btside);
            hCorrSideMassProj->SetLineColor(kBlack);
            
            TH1D *hCorrSideSideProj = hCorrSideSide[it][ia]->ProjectionX();
            hCorrSideSideProj->Scale(1.0/nPairSideSide);
            hCorrSideSideProj->Scale((betat*betaa)/(alphat*alphaa));
            //hCorrSideSideProj->Scale((bt*ba)/(btside*baside));
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
            hRatio->Draw();
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
    double mu1 = p[2];
    double mu2 = p[3];
    double sigma1 = p[4];
    double sigma2 = p[5];

    return c1*TMath::Exp(-(x[0]-mu1)*(x[0]-mu1)/(2*sigma1*sigma1)) + c2*TMath::Exp(-(x[0]-mu2)*(x[0]-mu2)/(2*sigma2*sigma2));
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
