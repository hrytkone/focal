void CompileEvent()
{
    if (!TClass::GetDict("AliJHMRCluster")) gROOT->ProcessLine(".L AliJHMRCluster.cxx++");
    if (!TClass::GetDict("AliJBaseTrack")) gROOT->ProcessLine(".L AliJBaseTrack.cxx++");
    if (!TClass::GetDict("AliJHMREvent")) gROOT->ProcessLine(".L AliJHMREvent.cxx++");
}
