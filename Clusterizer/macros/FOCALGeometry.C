AliFOCALGeometry * geometry;
AliRunLoader * fRunLoader;
AliFOCAL *fFOCAL;

void FOCALGeometry()
{
    gSystem->Load("/home/heimarry/alice/sw/ubuntu1804_x86-64/AliRoot/latest/lib/libpythia6.so");
    gSystem->Load("/home/heimarry/alice/sw/ubuntu1804_x86-64/AliRoot/latest/lib/libAliPythia6.so");
    gSystem->Load("libGeom");

    gGeoManager->Import("geometry.root","");

    TGeoVolume *topvol = gGeoManager->GetTopVolume();
    topvol->SetActiveDaughters();
    TGeoVolume *vol = topvol->GetNode(76)->GetVolume();

    //topvol->PrintNodes();

    topvol->Draw("ogl");

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
