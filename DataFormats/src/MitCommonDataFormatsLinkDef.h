// $Id: MitCommonDataFormatsLinkDef.h,v 1.6 2010/01/18 14:32:01 bendavid Exp $

#ifndef MITCOMMON_DATAFORMATS_LINKDEF_H
#define MITCOMMON_DATAFORMATS_LINKDEF_H

#include "MitCommon/DataFormats/interface/Hist1DRat.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitCommon/DataFormats/interface/TH3DAsymErr.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Vect3C.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
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
#pragma link C++ typedef mithep::UIntPair;
#pragma link C++ typedef mithep::UIntBounds;

#pragma link C++ class mithep::Hist1DRat+;
#pragma link C++ class mithep::TH2DAsymErr+;
#pragma link C++ enum mithep::TH2DAsymErr::EErrType;
#pragma link C++ class mithep::TH3DAsymErr+;
#pragma link C++ class mithep::Vect3+;
#pragma link C++ class mithep::Vect3C+;
#pragma link C++ class mithep::Vect4M+;
#endif
