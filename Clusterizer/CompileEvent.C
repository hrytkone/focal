void CompileEvent()
{
    //gROOT->LoadMacro("Event.cxx++g");
    if (!TClass::GetDict("AliJHMREvent")) {
        gROOT->ProcessLine(".L AliJHMREvent.cxx++");
    }
}
