#ifndef MITCOMMON_UTILS_LINKDEF_H
#define MITCOMMON_UTILS_LINKDEF_H

#include "MitCommon/Utils/interface/Utils.h"
#endif

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::Utils;
#endif
