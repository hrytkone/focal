void RunAcceptanceCheck(TString infile)
{
    gSystem->Load("../AliJBaseTrack_cxx.so");
    gSystem->Load("../AliJHMRCluster_cxx.so");
    gSystem->Load("../AliJHMREvent_cxx.so");
    gROOT->ProcessLine(Form(".x CheckFoCalAcceptance.C(\"%s\")", infile.Data()));
}
