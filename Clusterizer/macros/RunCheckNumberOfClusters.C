void RunCheckNumberOfClusters(TString infile)
{
    gSystem->Load("../src/AliJBaseTrack_cxx.so");
    gSystem->Load("../src/AliJHMRCluster_cxx.so");
    gSystem->Load("../src/AliJHMREvent_cxx.so");
    gROOT->ProcessLine(Form(".x CheckNumberOfClustersForGammas.C(\"%s\")", infile.Data()));
}
