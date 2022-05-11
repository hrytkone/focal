void RunAsymmetry(TString input)
{
    gSystem->Load("AliJBaseTrack_cxx.so");
    gSystem->Load("AliJHMRCluster_cxx.so");
    gSystem->Load("AliJHMREvent_cxx.so");
    gROOT->ProcessLine(Form(".x PlotAsymmetry.C(\"%s\")", input.Data()));
}
