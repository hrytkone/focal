void RunMassHistos(TString infile)
{
    gSystem->Load("../AliJBaseTrack_cxx.so");
    gSystem->Load("../AliJHMRCluster_cxx.so");
    gSystem->Load("../AliJHMREvent_cxx.so");
    gROOT->ProcessLine(Form(".x CreateMassHistos.C(\"%s\")", infile.Data()));
}
