AliFOCALGeometry * geometry;
AliRunLoader * fRunLoader;
AliFOCAL *fFOCAL;

void FOCALGeometry()
{
    gSystem->Load("/home/heimarry/alice/sw/ubuntu1804_x86-64/AliRoot/latest/lib/libpythia6.so");
    gSystem->Load("/home/heimarry/alice/sw/ubuntu1804_x86-64/AliRoot/latest/lib/libAliPythia6.so");
    gSystem->Load("libGeom");

    gGeoManager->Import("geometry.root","");

    TObjArray* volumes = gGeoManager->GetListOfVolumes();
    for (Int_t i=0; i<volumes->GetEntriesFast(); i++) {
        if ( !((TGeoVolume*)volumes->At(i))->IsAssembly() )
        ((TGeoVolume*)volumes->At(i))->SetVisibility(kTRUE);
    }


    //gGeoManager->SetVisLevel(5);
    //gGeoManager->GetTopVolume()->SetVisContainers(kTRUE);
    //gGeoManager->GetTopVolume()->Draw("ogl");

    TGeoVolume *topvol = gGeoManager->GetTopVolume();
    topvol->SetActiveDaughters();
    TGeoVolume *vol_focal = topvol->GetNode("FOCAL_1")->GetVolume();
    //TGeoVolume *vol_pipe = topvol->GetNode("RB26Pipe_1")->GetVolume();
    //vol_pipe->SetLineColor(kRed);

    gGeoManager->CloseGeometry();

    //topvol->PrintNodes();

    vol_focal->Draw("ogl");
    //vol_pipe->Draw("ogl");
    //gGeoManager->Edit();

//    geom->CloseGeometry();

    /**new TGeoManager("world", "the simplest geometry");

    TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
    TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
    TGeoVolume *top=gGeoManager->MakeBox("Top",med,10.,10.,10.);

    gGeoManager->SetTopVolume(top);
    gGeoManager->CloseGeometry();

    top->SetLineColor(kMagenta);
    gGeoManager->SetTopVisible(); // the TOP is invisible
    top->Draw();

    if (fRunLoader) fRunLoader->Delete();
    fRunLoader = AliRunLoader::Open("galice.root", "READ");

    if (!fRunLoader) {
        cout << "ERROR : RunLoader not found"  << endl;
    }

    if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
    gAlice = fRunLoader->GetAliRun();
    fFOCAL  = (AliFOCAL*)gAlice->GetDetector("FOCAL");**/

    //fFOCAL->Draw();
}
