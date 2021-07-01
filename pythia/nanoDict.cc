// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME nanoDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "include/AliJBaseTrack.h"
#include "include/AliJHMRHist.h"
#include "include/AliJHMRCorr.h"
#include "include/AliJHMRPythiaCatalyst.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *AliJBaseTrack_Dictionary();
   static void AliJBaseTrack_TClassManip(TClass*);
   static void *new_AliJBaseTrack(void *p = 0);
   static void *newArray_AliJBaseTrack(Long_t size, void *p);
   static void delete_AliJBaseTrack(void *p);
   static void deleteArray_AliJBaseTrack(void *p);
   static void destruct_AliJBaseTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJBaseTrack*)
   {
      ::AliJBaseTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJBaseTrack));
      static ::ROOT::TGenericClassInfo 
         instance("AliJBaseTrack", "include/AliJBaseTrack.h", 32,
                  typeid(::AliJBaseTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJBaseTrack_Dictionary, isa_proxy, 0,
                  sizeof(::AliJBaseTrack) );
      instance.SetNew(&new_AliJBaseTrack);
      instance.SetNewArray(&newArray_AliJBaseTrack);
      instance.SetDelete(&delete_AliJBaseTrack);
      instance.SetDeleteArray(&deleteArray_AliJBaseTrack);
      instance.SetDestructor(&destruct_AliJBaseTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJBaseTrack*)
   {
      return GenerateInitInstanceLocal((::AliJBaseTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJBaseTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJBaseTrack_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJBaseTrack*)0x0)->GetClass();
      AliJBaseTrack_TClassManip(theClass);
   return theClass;
   }

   static void AliJBaseTrack_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliJHMRHist_Dictionary();
   static void AliJHMRHist_TClassManip(TClass*);
   static void *new_AliJHMRHist(void *p = 0);
   static void *newArray_AliJHMRHist(Long_t size, void *p);
   static void delete_AliJHMRHist(void *p);
   static void deleteArray_AliJHMRHist(void *p);
   static void destruct_AliJHMRHist(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJHMRHist*)
   {
      ::AliJHMRHist *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJHMRHist));
      static ::ROOT::TGenericClassInfo 
         instance("AliJHMRHist", "include/AliJHMRHist.h", 18,
                  typeid(::AliJHMRHist), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJHMRHist_Dictionary, isa_proxy, 0,
                  sizeof(::AliJHMRHist) );
      instance.SetNew(&new_AliJHMRHist);
      instance.SetNewArray(&newArray_AliJHMRHist);
      instance.SetDelete(&delete_AliJHMRHist);
      instance.SetDeleteArray(&deleteArray_AliJHMRHist);
      instance.SetDestructor(&destruct_AliJHMRHist);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJHMRHist*)
   {
      return GenerateInitInstanceLocal((::AliJHMRHist*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJHMRHist*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJHMRHist_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJHMRHist*)0x0)->GetClass();
      AliJHMRHist_TClassManip(theClass);
   return theClass;
   }

   static void AliJHMRHist_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliJHMRCorr_Dictionary();
   static void AliJHMRCorr_TClassManip(TClass*);
   static void *new_AliJHMRCorr(void *p = 0);
   static void *newArray_AliJHMRCorr(Long_t size, void *p);
   static void delete_AliJHMRCorr(void *p);
   static void deleteArray_AliJHMRCorr(void *p);
   static void destruct_AliJHMRCorr(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJHMRCorr*)
   {
      ::AliJHMRCorr *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJHMRCorr));
      static ::ROOT::TGenericClassInfo 
         instance("AliJHMRCorr", "include/AliJHMRCorr.h", 26,
                  typeid(::AliJHMRCorr), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJHMRCorr_Dictionary, isa_proxy, 0,
                  sizeof(::AliJHMRCorr) );
      instance.SetNew(&new_AliJHMRCorr);
      instance.SetNewArray(&newArray_AliJHMRCorr);
      instance.SetDelete(&delete_AliJHMRCorr);
      instance.SetDeleteArray(&deleteArray_AliJHMRCorr);
      instance.SetDestructor(&destruct_AliJHMRCorr);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJHMRCorr*)
   {
      return GenerateInitInstanceLocal((::AliJHMRCorr*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJHMRCorr*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJHMRCorr_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJHMRCorr*)0x0)->GetClass();
      AliJHMRCorr_TClassManip(theClass);
   return theClass;
   }

   static void AliJHMRCorr_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliJHMRPythiaCatalyst_Dictionary();
   static void AliJHMRPythiaCatalyst_TClassManip(TClass*);
   static void delete_AliJHMRPythiaCatalyst(void *p);
   static void deleteArray_AliJHMRPythiaCatalyst(void *p);
   static void destruct_AliJHMRPythiaCatalyst(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJHMRPythiaCatalyst*)
   {
      ::AliJHMRPythiaCatalyst *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJHMRPythiaCatalyst));
      static ::ROOT::TGenericClassInfo 
         instance("AliJHMRPythiaCatalyst", "include/AliJHMRPythiaCatalyst.h", 29,
                  typeid(::AliJHMRPythiaCatalyst), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJHMRPythiaCatalyst_Dictionary, isa_proxy, 0,
                  sizeof(::AliJHMRPythiaCatalyst) );
      instance.SetDelete(&delete_AliJHMRPythiaCatalyst);
      instance.SetDeleteArray(&deleteArray_AliJHMRPythiaCatalyst);
      instance.SetDestructor(&destruct_AliJHMRPythiaCatalyst);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJHMRPythiaCatalyst*)
   {
      return GenerateInitInstanceLocal((::AliJHMRPythiaCatalyst*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJHMRPythiaCatalyst*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJHMRPythiaCatalyst_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJHMRPythiaCatalyst*)0x0)->GetClass();
      AliJHMRPythiaCatalyst_TClassManip(theClass);
   return theClass;
   }

   static void AliJHMRPythiaCatalyst_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJBaseTrack(void *p) {
      return  p ? new(p) ::AliJBaseTrack : new ::AliJBaseTrack;
   }
   static void *newArray_AliJBaseTrack(Long_t nElements, void *p) {
      return p ? new(p) ::AliJBaseTrack[nElements] : new ::AliJBaseTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJBaseTrack(void *p) {
      delete ((::AliJBaseTrack*)p);
   }
   static void deleteArray_AliJBaseTrack(void *p) {
      delete [] ((::AliJBaseTrack*)p);
   }
   static void destruct_AliJBaseTrack(void *p) {
      typedef ::AliJBaseTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJBaseTrack

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJHMRHist(void *p) {
      return  p ? new(p) ::AliJHMRHist : new ::AliJHMRHist;
   }
   static void *newArray_AliJHMRHist(Long_t nElements, void *p) {
      return p ? new(p) ::AliJHMRHist[nElements] : new ::AliJHMRHist[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJHMRHist(void *p) {
      delete ((::AliJHMRHist*)p);
   }
   static void deleteArray_AliJHMRHist(void *p) {
      delete [] ((::AliJHMRHist*)p);
   }
   static void destruct_AliJHMRHist(void *p) {
      typedef ::AliJHMRHist current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJHMRHist

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJHMRCorr(void *p) {
      return  p ? new(p) ::AliJHMRCorr : new ::AliJHMRCorr;
   }
   static void *newArray_AliJHMRCorr(Long_t nElements, void *p) {
      return p ? new(p) ::AliJHMRCorr[nElements] : new ::AliJHMRCorr[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJHMRCorr(void *p) {
      delete ((::AliJHMRCorr*)p);
   }
   static void deleteArray_AliJHMRCorr(void *p) {
      delete [] ((::AliJHMRCorr*)p);
   }
   static void destruct_AliJHMRCorr(void *p) {
      typedef ::AliJHMRCorr current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJHMRCorr

namespace ROOT {
   // Wrapper around operator delete
   static void delete_AliJHMRPythiaCatalyst(void *p) {
      delete ((::AliJHMRPythiaCatalyst*)p);
   }
   static void deleteArray_AliJHMRPythiaCatalyst(void *p) {
      delete [] ((::AliJHMRPythiaCatalyst*)p);
   }
   static void destruct_AliJHMRPythiaCatalyst(void *p) {
      typedef ::AliJHMRPythiaCatalyst current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJHMRPythiaCatalyst

namespace {
  void TriggerDictionaryInitialization_nanoDict_Impl() {
    static const char* headers[] = {
"include/AliJBaseTrack.h",
"include/AliJHMRHist.h",
"include/AliJHMRCorr.h",
"include/AliJHMRPythiaCatalyst.h",
0
    };
    static const char* includePaths[] = {
"/home/alidock/.sw/slc7_x86-64/pythia/latest/include",
"/home/alidock/.sw/slc7_x86-64/ROOT/v6-20-08-alice1-73/include/",
"/mnt/Desktop/focal/pythia/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "nanoDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$include/AliJBaseTrack.h")))  AliJBaseTrack;
class __attribute__((annotate("$clingAutoload$include/AliJHMRHist.h")))  AliJHMRHist;
class __attribute__((annotate("$clingAutoload$include/AliJHMRCorr.h")))  AliJHMRCorr;
class __attribute__((annotate("$clingAutoload$include/AliJHMRPythiaCatalyst.h")))  AliJHMRPythiaCatalyst;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "nanoDict dictionary payload"

#ifndef JTKT
  #define JTKT 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "include/AliJBaseTrack.h"
#include "include/AliJHMRHist.h"
#include "include/AliJHMRCorr.h"
#include "include/AliJHMRPythiaCatalyst.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"AliJBaseTrack", payloadCode, "@",
"AliJHMRCorr", payloadCode, "@",
"AliJHMRHist", payloadCode, "@",
"AliJHMRPythiaCatalyst", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("nanoDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_nanoDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_nanoDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_nanoDict() {
  TriggerDictionaryInitialization_nanoDict_Impl();
}
