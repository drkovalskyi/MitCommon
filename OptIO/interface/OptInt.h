//--------------------------------------------------------------------------------------------------
// $Id: OptInt.h,v 1.1 2009/02/24 11:56:42 loizides Exp $
//
// OptInt
//
// Interface to set some of the optimized IO parameters.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_OPTIO_OPTINT_H
#define MITCOMMON_OPTIO_OPTINT_H

#include <TObject.h>
 
namespace mithep 
{
  class OptInt
  {
    public:
      OptInt() {}

      static void SetAlgoFractions(Double_t lzo, Double_t gz, Double_t bz, Double_t lzma);
      static void SetBzipFraction(Double_t f);
      static void SetGzipFraction(Double_t f);
      static void SetLzmaFraction(Double_t f);
      static void SetLzoFraction(Double_t f);
      static void SetStaticMalloc(Bool_t b);
      static void SetVerbose(Int_t vl);
      static void SetZipMode(Int_t zm);

    ClassDef(OptInt, 0) // Interface to OptIO parameters
  };
}
#endif
