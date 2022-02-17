void RunClusterizer(TString indir, Int_t njobs)
{
    gSystem->Load("AliJHMREvent_cxx.so");
    gROOT->ProcessLine(Form(".x Clusterizer.C(\"%s\",%d)", indir.Data(), njobs));
}
