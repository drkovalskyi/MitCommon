// $Id: MitCommonDataFormatsLinkDef.h,v 1.3 2009/04/29 15:14:06 loizides Exp $

#ifndef MITCOMMON_MATHTOOLS_LINKDEF_H
#define MITCOMMON_MATHTOOLS_LINKDEF_H

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/MathTools/interface/TPFuncs.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::MathUtils-;
#pragma link C++ class mithep::TPFuncs-;
#endif
