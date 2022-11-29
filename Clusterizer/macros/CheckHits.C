AliFOCALGeometry * fGeom;
AliRunLoader * fRunLoader;
AliFOCAL *fFOCAL;
AliFOCALLoader *fFOCALLoader;
TDirectory * actualDir;

TFile *outfile;
TBranch * branchHits;

const int nlayer = 20;
const int nseg = 6;
TH2D *hHitEnergyMap[nlayer];
TH2D *hHitEnergyMapSeg[nseg];

void InitOutput(TString outputname);
void LoadClusterizerHit(TString inputfile);

void CheckHits(TString simFolder, int ifolder, TString outputname)
{
    gSystem->Load("/home/heimarry/alice/sw/ubuntu1804_x86-64/AliRoot/latest/lib/libpythia6.so");
    gSystem->Load("/home/heimarry/alice/sw/ubuntu1804_x86-64/AliRoot/latest/lib/libAliPythia6.so");

    fGeom = AliFOCALGeometry::GetInstance("../geometry.txt");

    InitOutput(outputname);

    int nevtot = 0;

    TString inputFile = Form("%s/%03d/%s", simFolder.Data(), ifolder, "galice.root");
    LoadClusterizerHit(inputFile);

    int nev = fRunLoader->GetNumberOfEvents();
    //nev = 1;
    cout << "\tN EVENTS : " << nev << endl;
    for (int ievt = 0; ievt < nev; ievt++) {
        fRunLoader->GetEvent(ievt);
        TTree* treeH = fFOCALLoader->TreeH();
        branchHits = treeH->GetBranch("FOCAL");

        TClonesArray* hits = 0;
        branchHits->SetAddress(&hits);

        for (Int_t iTrack=0; iTrack<branchHits->GetEntries(); iTrack++) {
            branchHits->GetEntry(iTrack); // load iTrack

            int nfocal = hits->GetEntries();
            // loop over hits
            for (Int_t ifocal = 0; ifocal < nfocal; ifocal++) {
                AliFOCALhit * fFOCALHit = (AliFOCALhit*) hits->UncheckedAt(ifocal);
                Float_t x, y, z, e;
                x = (Float_t) fFOCALHit->X();
                y = (Float_t) fFOCALHit->Y();
                z = (Float_t) fFOCALHit->Z();
                e = (Float_t) fFOCALHit->GetEnergy();

                Int_t col,row,layer,seg;
                fGeom->GetVirtualInfo(x,y,z,col,row,layer,seg);
                //if (seg==1) cout << "layer for segement 1 : " << layer << endl;

                if (layer > -1 && layer < 20)
                    hHitEnergyMap[layer]->Fill(x,y,e);

                if (seg > -1 && seg < 6)
                    hHitEnergyMapSeg[seg]->Fill(x,y,e);

            } // Hits Loop ended

        }
    } // end of event loop
    outfile->cd();
    outfile->Write("", TObject::kOverwrite);
}

//******************************************************************************
//******************************************************************************
void InitOutput(TString outputname)
{
    outfile = TFile::Open(Form("%s.root", outputname.Data()), "RECREATE");
    for (int i=0; i<nlayer; i++)
        hHitEnergyMap[i] = new TH2D(Form("hHitEnergyMap_L%d", i), Form("Hit energy map layer %d", i), 500, -50., 50., 500, -50., 50.);

    for (int i=0; i<nseg; i++)
        hHitEnergyMapSeg[i] = new TH2D(Form("hHitEnergyMapSeg_L%d", i), Form("Hit energy map segment %d", i), 500, -50., 50., 500, -50., 50.);
}

void LoadClusterizerHit(TString inputfile)
{
    if (fRunLoader) fRunLoader->Delete();
    fRunLoader = AliRunLoader::Open(inputfile.Data(), "READ");

    if (!fRunLoader) {
        cout << "ERROR : RunLoader not found for " << inputfile.Data() << endl;
    }

    if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
    if (!fRunLoader->TreeE()) fRunLoader->LoadHeader();
    if (!fRunLoader->TreeK()) fRunLoader->LoadKinematics();

    gAlice = fRunLoader->GetAliRun();

    fFOCAL  = (AliFOCAL*)gAlice->GetDetector("FOCAL");
    fFOCALLoader = dynamic_cast<AliFOCALLoader*>(fRunLoader->GetLoader("FOCALLoader"));
    fFOCALLoader->LoadHits("READ");
}
