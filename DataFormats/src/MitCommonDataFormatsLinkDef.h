// $Id: MitAnaDataTreeLinkDef.h,v 1.38 2008/09/19 11:58:03 bendavid Exp $

#ifndef MITCOMMON_DATAFORMATS_LINKDEF_H
#define MITCOMMON_DATAFORMATS_LINKDEF_H

#include "MitCommon/DataFormats/interface/Types.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ typedef mithep::FourVector;
#pragma link C++ typedef mithep::FourVectorM;
#pragma link C++ typedef mithep::FourVectorE;
#pragma link C++ typedef mithep::ThreeVector;
#pragma link C++ typedef mithep::ThreeVectorC;
#pragma link C++ typedef mithep::ThreeSymMatrix;
#pragma link C++ typedef mithep::SevenSymMatrix;
#pragma link C++ typedef mithep::ThreeMatrix;
#pragma link C++ typedef mithep::SevenMatrix;
#pragma link C++ typedef mithep::FourVector32;
#pragma link C++ typedef mithep::FourVectorM32;
#pragma link C++ typedef mithep::FourVectorE32;
#pragma link C++ typedef mithep::ThreeVector32;
#pragma link C++ typedef mithep::ThreeVectorC32;
#pragma link C++ typedef mithep::ThreeSymMatrix32;
#pragma link C++ typedef mithep::SevenSymMatrix32;
#pragma link C++ typedef mithep::ThreeMatrix32;
#pragma link C++ typedef mithep::SevenMatrix32;

#endif