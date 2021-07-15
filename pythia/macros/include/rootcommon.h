void hset(TH1& hid, TString xtit="", TString ytit="",
        double titoffx = 1.1, double titoffy = 1.1,
        double titsizex = 0.06, double titsizey = 0.06,
        double labeloffx = 0.01, double labeloffy = 0.001,
        double labelsizex = 0.05, double labelsizey = 0.05,
        int divx = 505, int divy=505)
{
    hid.GetXaxis()->CenterTitle(1);
    hid.GetYaxis()->CenterTitle(1);

    hid.GetXaxis()->SetTitleOffset(titoffx);
    hid.GetYaxis()->SetTitleOffset(titoffy);

    hid.GetXaxis()->SetTitleSize(titsizex);
    hid.GetYaxis()->SetTitleSize(titsizey);

    hid.GetXaxis()->SetLabelOffset(labeloffx);
    hid.GetYaxis()->SetLabelOffset(labeloffy);

    hid.GetXaxis()->SetLabelSize(labelsizex);
    hid.GetYaxis()->SetLabelSize(labelsizey);

    hid.GetXaxis()->SetNdivisions(divx);
    hid.GetYaxis()->SetNdivisions(divy);

    hid.GetXaxis()->SetTitle(xtit);
    hid.GetYaxis()->SetTitle(ytit);

    hid.GetXaxis()->SetLabelFont(42);
    hid.GetYaxis()->SetLabelFont(42);
    hid.GetXaxis()->SetTitleFont(42);
    hid.GetYaxis()->SetTitleFont(42);
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
            char name[16];
            sprintf(name,"pad_%i_%i",i,j);
            TPad *pad = (TPad*) gROOT->FindObject(name);
            if (pad) delete pad;
            pad = new TPad(name,"",hposl,vposd,hposr,vposu);
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

TGraphErrors *get_graph_total_error( TGraphErrors  *stat, TGraphErrors *syst){
    TGraphErrors * gr_total_error = new TGraphErrors( stat->GetN() );
    //TGraph graph_error_syst( syst->GetN(),  syst->GetX(),syst->GetEY() );
    for( int i=0; i< stat->GetN(); i++ ){
        double x = stat->GetX()[i];
        double y = stat->GetY()[i];
        double ey2 = stat->GetErrorY(i);
        double ey1 = syst->GetErrorY(i);
        //double ey1 = graph_error_syst.Eval(x);
        gr_total_error->SetPoint( i,  x, y);
        gr_total_error->SetPointError( i,  0, TMath::Sqrt( ey1*ey1+ey2*ey2));
    }
    return gr_total_error;
}

TGraphErrors *get_ratio( TGraphErrors  *l, TGraphErrors *r){
    TGraphErrors * gr_ratio = new TGraphErrors( r->GetN() );
    TGraph gel( l->GetN(),  l->GetX(),l->GetEY() );
    for( int i=0; i< r->GetN(); i++ ){
        double x = r->GetX()[i];
        double y2 = r->GetY()[i];
        //double ey2 = r->GetEY()[i];
        double ey2 = r->GetErrorY(i);
        double y1 = l->Eval(x);
        double ey1 = gel.Eval(x);


        double ratio = y1 / y2;
        gr_ratio->SetPoint( i,  x, ratio);
        gr_ratio->SetPointError( i,  0, ratio*TMath::Sqrt( ey1*ey1/y1/y1+ey2*ey2/y2/y2));
    }
    gr_ratio->SetMarkerStyle(l->GetMarkerStyle());
    gr_ratio->SetMarkerColor(l->GetMarkerColor());
    gr_ratio->SetLineColor(l->GetLineColor());

    return gr_ratio;
}


TGraphAsymmErrors *GetDataOverTheory(TGraphAsymmErrors *gr, TF1 *ftheory ){
    double x[300], y[300], exl[300], exh[300], eyl[300], eyh[300];
    int NC =  gr->GetN();
    for(int ii=0;ii<NC;ii++){
        gr->GetPoint(ii,x[ii],y[ii]);
        exl[ii] = gr->GetErrorXlow(ii);
        exh[ii] = gr->GetErrorXhigh(ii);
        eyl[ii] = gr->GetErrorYlow(ii);
        eyh[ii] = gr->GetErrorYhigh(ii);
    }
    for(int ii=0;ii<NC;ii++){
        //y[ii]   = ( ftheory->Eval(x[ii]) - y[ii] )/y[ii];
        //eyl[ii] = ftheory->Eval(x[ii])*eyl[ii] / y[ii] / y[ii];
        //eyh[ii] = ftheory->Eval(x[ii])*eyh[ii] / y[ii] / y[ii];
        y[ii]   = ( ftheory->Eval(x[ii]) - y[ii] ) / ftheory->Eval(x[ii]);
        eyl[ii] = eyl[ii]/ftheory->Eval(x[ii]);
        eyh[ii] = eyh[ii]/ftheory->Eval(x[ii]);
    }
    return new TGraphAsymmErrors(NC, x, y, exl, exh, eyl, eyh);
}

TGraphAsymmErrors *GetDataOverTheory(TGraphAsymmErrors *grData, TGraphAsymmErrors *grTheo ){
    double x[300], y[300], exl[300], exh[300], eyl[300], eyh[300];
    int NC =  grData->GetN();
    for(int ii=0;ii<NC;ii++){
        grData->GetPoint(ii,x[ii],y[ii]);
        exl[ii] = grData->GetErrorXlow(ii);
        exh[ii] = grData->GetErrorXhigh(ii);
        eyl[ii] = grData->GetErrorYlow(ii);
        eyh[ii] = grData->GetErrorYhigh(ii);
    }
    for(int ii=0;ii<NC;ii++){
        y[ii]   = ( grTheo->Eval(x[ii]) - y[ii] ) / grTheo->Eval(x[ii]);
        eyl[ii] = eyl[ii]/grTheo->Eval(x[ii]);
        eyh[ii] = eyh[ii]/grTheo->Eval(x[ii]);
    }
    return new TGraphAsymmErrors(NC, x, y, exl, exh, eyl, eyh);
}

TGraphErrors* grrScale( TGraphErrors *gr, double sc=1){
    double xin, yin, exin, eyin;
    double x[100], y[100], ex[100]={0}, ey[100];
    int NC =  gr->GetN();
    for(int ii=0;ii<NC;ii++){
        gr->GetPoint(ii, xin, yin); exin = gr->GetErrorX(ii); eyin = gr->GetErrorY(ii);
        x[ii] = xin;
        y[ii] = yin * sc;
        ex[ii]= exin;
        ey[ii]= eyin * sc;
    }
    gr = new TGraphErrors(NC, x, y, ex, ey);
    return gr;
}

void scaleThisGrr( TGraphErrors &gr, double sc=1){
    int NC =  gr.GetN();
    double xin, yin, exin, eyin;
    for(int ii=0;ii<NC;ii++){
        gr.GetPoint(ii, xin, yin);    exin = gr.GetErrorX(ii); eyin = gr.GetErrorY(ii);
        gr.SetPoint(ii, xin, yin*sc); gr.SetPointError(ii, exin, eyin*sc);
    }
}

void scaleThisGrr( TGraphAsymmErrors &gr, double sc=1){
    int NC =  gr.GetN();
    double xin, yin, exinl, exinh, eyinl, eyinh;
    for(int ii=0;ii<NC;ii++){
        gr.GetPoint(ii, xin, yin);    
        exinl = gr.GetErrorXlow(ii);  exinh = gr.GetErrorXhigh(ii); 
        eyinl = gr.GetErrorYlow(ii);  eyinh = gr.GetErrorYhigh(ii); 
        gr.SetPoint(ii, xin, yin*sc); gr.SetPointError(ii, exinl, exinh, eyinl*sc, eyinh*sc);
    }
}

//This calculates the chi^2/Ndof for the experiment vs theory. The uncertainty used is sigma^2 = sigma_{exp_stat}^2 + sigma_{exp_syst}^2+sigma_{theory}^2
double chi_squared(TGraphErrors *expStat, TGraphErrors *expSyst, TGraphErrors *theory) {
    double sum = 0;
    TGraph gel( theory->GetN(),  theory->GetX(),theory->GetEY() );
    for(int i= 0; i<expStat->GetN(); i++) {
        double x = expStat->GetX()[i];
        double eyStat = expStat->GetErrorY(i);
        double eySyst = expSyst->GetErrorY(i);
        double eyTheory = gel.Eval(x);
        double diff =(expStat->GetY()[i] - theory->Eval(x) ); 
        sum += diff*diff /( eyStat*eyStat + eySyst*eySyst +eyTheory*eyTheory );
    }
    return (sum/expStat->GetN());
}

void SetTGraphXError(TGraphErrors *gr,double xerr){
    int NC =  gr->GetN();
    for(int ii=0;ii<NC;ii++){
        gr->SetPointError(ii,xerr,gr->GetErrorY(ii));
    }
}
