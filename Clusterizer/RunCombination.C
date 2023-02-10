void RunCombination(TString simfolder, TString clusterfolder, TString outputdir, Int_t startFolder, Int_t endFolder)
{
    gSystem->Load("src/AliJBaseTrack_cxx.so");
    gSystem->Load("src/AliJHMRCluster_cxx.so");
    gSystem->Load("src/AliJHMREvent_cxx.so");
    gROOT->ProcessLine(Form(".x CombineSimulationOutput.C(\"%s\",\"%s\",\"%s\",%d,%d)", simfolder.Data(), clusterfolder.Data(), outputdir.Data(), startFolder, endFolder));
}
