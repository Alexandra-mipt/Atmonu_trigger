// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME AtmoNuTrigger_dict

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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "classes.h"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_libAtmoNuTrigger_dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libAtmoNuTrigger_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libAtmoNuTrigger_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef _POSIX_SOURCE
  #define _POSIX_SOURCE 1
#endif
#ifndef _SVID_SOURCE
  #define _SVID_SOURCE 1
#endif
#ifndef _BSD_SOURCE
  #define _BSD_SOURCE 1
#endif
#ifndef _POSIX_C_SOURCE
  #define _POSIX_C_SOURCE 2
#endif
#ifndef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 1
#endif
#ifndef UNIX
  #define UNIX 1
#endif
#ifndef LINUX
  #define LINUX 1
#endif
#ifndef __UNIX__
  #define __UNIX__ 1
#endif
#ifndef __LINUX__
  #define __LINUX__ 1
#endif
#ifndef DATAREP_LITTLE_IEEE
  #define DATAREP_LITTLE_IEEE 1
#endif
#ifndef DATAREP_LITTLE_ENDIAN
  #define DATAREP_LITTLE_ENDIAN 1
#endif
#ifndef DEFECT_NO_IOSTREAM_NAMESPACES
  #define DEFECT_NO_IOSTREAM_NAMESPACES 1
#endif
#ifndef DEFECT_NO_JZEXT
  #define DEFECT_NO_JZEXT 1
#endif
#ifndef DEFECT_NO_INTHEX
  #define DEFECT_NO_INTHEX 1
#endif
#ifndef DEFECT_NO_INTHOLLERITH
  #define DEFECT_NO_INTHOLLERITH 1
#endif
#ifndef DEFECT_NO_READONLY
  #define DEFECT_NO_READONLY 1
#endif
#ifndef DEFECT_NO_DIRECT_FIXED
  #define DEFECT_NO_DIRECT_FIXED 1
#endif
#ifndef DEFECT_NO_STRUCTURE
  #define DEFECT_NO_STRUCTURE 1
#endif
#ifndef USE_ROOT
  #define USE_ROOT 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/nova/alexandra/rel_development/DAQDataFormats/cxx/include/TriggerDefines.h"



#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"daqdataformats::trigID", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libAtmoNuTrigger_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libAtmoNuTrigger_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libAtmoNuTrigger_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libAtmoNuTrigger_dict() {
  TriggerDictionaryInitialization_libAtmoNuTrigger_dict_Impl();
}
