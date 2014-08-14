//--------------------------------------------------------------------------------------------------
// $Id: OptInt.h,v 1.5 2009/03/04 07:24:49 loizides Exp $
//
// OptInt
//
// Interface to set some of the optimized IO parameters. Use Compress or Decompress 
// to compress or decompress data.
//
// The following zip modes are supported:
//   Set ZipMode to be used in R__zip. Supported values range from 1 to 5, where
//    1 == gzip  (http://www.gzip.org/, standard in ROOT)
//    2 == bzip2 (http://www.bzip.org/)
//    3 == lzo   (http://www.oberhumer.com/opensource/lzo/)
//    4 == rle   (http://bcl.comli.eu/)
//    5 == lzma  (http://www.7-zip.org/sdk.html)
//    99 == combination of 1 to 5.
// 
// Additionally you can specify the compression fraction that must be reached to accept result of 
// the "heavy" compression algorithms. Negative values turn off usage of the specific compression 
// algorithm in ZipMode=99.
//
// In order to make use of all of this, make sure you preload MitCommonOptIO.so. This can be 
// achieved by setting LD_PRELOAD to $CMSSW_BASE/lib/$SCRAM_ARCH/libMitCommonOptIO.so.
// 
// In additon, you can use .rootrc, for example like
//  Root.ZipMode: 5
//  Root.OptIO.Verbose: 0
//  Root.OptIO.SMalloc: 0
//  Root.OptIO.LzoFraction:  1
//  Root.OptIO.GzipFraction: 1
//  Root.OptIO.BzipFraction: 1
//  Root.OptIO.LzmaFraction: 1
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
      virtual ~OptInt() {}

      static Int_t  Compress(Int_t srcsize, char *src, Int_t tgtsize, char *tgt, 
                             Int_t cl, Int_t method=-1);
      static Int_t  DeCompress(Int_t srcsize, char *src, Int_t tgtsize, char *tgt); 
      static Bool_t IsActivated();
      static void   SetAlgoFractions(Double_t lzo, Double_t gz, Double_t bz, Double_t lzma);
      static void   SetBzipFraction(Double_t f);
      static void   SetGzipFraction(Double_t f);
      static void   SetLzmaFraction(Double_t f);
      static void   SetLzoFraction(Double_t f);
      static void   SetStaticMalloc(Bool_t b);
      static void   SetVerbose(Int_t vl);
      static void   SetZipMode(Int_t zm);

    ClassDef(OptInt, 0) // Interface to OptIO parameters
  };
}
#endif
