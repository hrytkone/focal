void CompileEvent()
{
    if (!TClass::GetDict("AliJHMRCluster")) gROOT->ProcessLine(".L src/AliJHMRCluster.cxx++");
    if (!TClass::GetDict("AliJBaseTrack")) gROOT->ProcessLine(".L src/AliJBaseTrack.cxx++");
    if (!TClass::GetDict("AliJHMREvent")) gROOT->ProcessLine(".L src/AliJHMREvent.cxx++");
}
